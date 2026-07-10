#include <vector>
#include <string>
#include <webp/encode.h>
#include <webp/decode.h>
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <thread>
#include <mutex>
#include <queue>
#include <condition_variable>
#include <atomic>
#include <filesystem>
#include <chrono>
#include <memory>
#include <functional>

namespace fs = std::filesystem;
struct WebpEncodeParams {
    int lossless;   // 1 = 无损，0 = 有损
    int quality;    // 0-100
    int method;     // 0-6
};

// ==================== 线程池 ====================
class ThreadPool {
public:
    ThreadPool(size_t threads) : stop(false) {
        for (size_t i = 0; i < threads; ++i)
            workers.emplace_back([this] {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            });
    }

    template<class F>
    void enqueue(F&& f) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread& worker : workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

// ==================== 核心：单块 WebP 编码（块级任务单元） ====================
// 供线程池调用，matrix 通过 shared_ptr 只读共享，零拷贝
void EncodeSingleWebpBlock(
    const std::vector<std::vector<uint8_t>>& matrix,
    const std::string& output_dir,
    int read_length,
    size_t block_idx,
    size_t start_seq,
    size_t end_seq,
    int image_width,
    int sequences_per_row,
    WebpEncodeParams webp_params)
{
    const size_t seqs_in_block = end_seq - start_seq;
    if (seqs_in_block == 0) return;

    int block_height = static_cast<int>(
        (seqs_in_block + sequences_per_row - 1) / sequences_per_row);

    // ARGB 缓冲区（4字节/像素）
    const int stride = image_width * 4;
    std::vector<uint8_t> argb_data(
        static_cast<size_t>(stride) * block_height, 0);

    // 填充 ARGB 数据
    for (size_t seq_offset = 0; seq_offset < seqs_in_block; ++seq_offset) {
        const size_t global_seq_idx = start_seq + seq_offset;
        if (global_seq_idx >= matrix.size()) break;

        const size_t row_in_img = seq_offset / sequences_per_row;
        const size_t col_in_row = seq_offset % sequences_per_row;
        const int    x_start    = static_cast<int>(col_in_row) * read_length;

        for (int pixel = 0; pixel < read_length; ++pixel) {
            const int x        = x_start + pixel;
            const int base_idx = static_cast<int>(row_in_img) * stride + x * 4;
            if (base_idx + 3 >= static_cast<int>(argb_data.size())) break;

            uint8_t v = matrix[global_seq_idx][pixel];
            argb_data[base_idx]     = v;
            argb_data[base_idx + 1] = v;
            argb_data[base_idx + 2] = v;
            argb_data[base_idx + 3] = 255; // Alpha 不透明
        }
    }

    // WebP 无损编码配置
    WebPConfig config;

    float q = static_cast<float>(std::clamp(webp_params.quality, 0, 100));

    if (!WebPConfigPreset(&config, WEBP_PRESET_DEFAULT, q)) {
        throw std::runtime_error("块" + std::to_string(block_idx) + " WebP配置初始化失败");
    }

    config.lossless = webp_params.lossless ? 1 : 0;
    config.method   = std::clamp(webp_params.method, 0, 6);
    config.quality  = q;

    // 对无损模式，强制 quality=100 更直观。
    // 如果你想让 lossless quality 作为压缩努力参数，也可以不强制。
    if (config.lossless) {
        config.quality = 100.0f;
    }

    if (!WebPValidateConfig(&config)) {
        throw std::runtime_error("块" + std::to_string(block_idx) + " WebP配置非法");
    }

    // 初始化图片结构
    WebPPicture pic;
    if (!WebPPictureInit(&pic))
        throw std::runtime_error("块" + std::to_string(block_idx) + " WebPPicture初始化失败");
    pic.width    = image_width;
    pic.height   = block_height;
    pic.use_argb = 1;

    if (!WebPPictureImportRGBA(&pic, argb_data.data(), stride)) {
        std::string err = "块" + std::to_string(block_idx) +
                          " RGBA导入失败，错误码：" + std::to_string(pic.error_code);
        WebPPictureFree(&pic);
        throw std::runtime_error(err);
    }

    // 内存写入器
    WebPMemoryWriter writer;
    WebPMemoryWriterInit(&writer);
    pic.writer     = WebPMemoryWrite;
    pic.custom_ptr = &writer;

    // 编码
    if (!WebPEncode(&config, &pic)) {
        std::string err = "块" + std::to_string(block_idx) +
                          " 无损编码失败，错误码：" + std::to_string(pic.error_code);
        WebPMemoryWriterClear(&writer);
        WebPPictureFree(&pic);
        throw std::runtime_error(err);
    }

    // 写文件
    const std::string filename = output_dir + "/lossless_block_" +
                                 std::to_string(block_idx) + ".webp";
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs || !ofs.write(reinterpret_cast<char*>(writer.mem), writer.size)) {
        WebPMemoryWriterClear(&writer);
        WebPPictureFree(&pic);
        throw std::runtime_error("块" + std::to_string(block_idx) + " 文件写入失败：" + filename);
    }

    WebPMemoryWriterClear(&writer);
    WebPPictureFree(&pic);
}

// ==================== 读取质量分块（不变） ====================
std::vector<std::string> ReadQualityScoreBlocks_40(int id,
                                                   char* quality_score_dir,
                                                   int read_length)
{
    std::string filename = std::string(quality_score_dir) +
                           "/quality_score." + std::to_string(id) + ".bin";

    std::ifstream in_file(filename, std::ios::binary);
    if (!in_file)
        throw std::runtime_error("无法打开文件: " + filename);

    size_t block_count = 0;
    in_file.read(reinterpret_cast<char*>(&block_count), sizeof(size_t));

    std::vector<std::string> blocks;
    blocks.reserve(block_count);

    for (size_t i = 0; i < block_count; ++i) {
        size_t code_len = 0;
        in_file.read(reinterpret_cast<char*>(&code_len), sizeof(size_t));

        std::string str(code_len, '\0');
        in_file.read(&str[0], code_len);

        if (code_len != static_cast<size_t>(read_length)) {
            std::cerr << "警告: 文件" << id << "第" << i
                      << "块长度错误，预期" << read_length
                      << "，实际" << code_len << std::endl;
        } else {
            blocks.emplace_back(std::move(str));
        }
    }

    if (in_file.fail() && !in_file.eof())
        throw std::runtime_error("读取文件失败: " + filename);
    in_file.close();

    return blocks;
}
void webp_main_40(char* quality_score_dir,
                  char* webp_dir,
                  int thread_n,
                  int total,
                  int read_length,
                  int webp_lossless,
                  int webp_quality,
                  int webp_method,
                  int webp_width,
                  int webp_height)
{
    try {
        const int TOTAL_FILES = total + 1;
        const int MAX_INFLIGHT_FILES = 50;  // ← 最多同时在队列里的文件数
        std::cout << "Webp压缩总共：" << TOTAL_FILES << std::endl;

        const int requested_width = std::clamp(webp_width, 1, WEBP_MAX_DIMENSION);
        const int requested_height = std::clamp(webp_height, 1, WEBP_MAX_DIMENSION);
        const int sequences_per_row   = requested_width / read_length;
        if (sequences_per_row == 0)
            throw std::runtime_error("read_length超过WebP最大宽度限制");
        const int    BASE_IMAGE_WIDTH    = sequences_per_row * read_length;
        const int    MAX_BLOCK_HEIGHT    = requested_height;
        const size_t max_seqs_per_block  =
            static_cast<size_t>(sequences_per_row) * MAX_BLOCK_HEIGHT;

        std::cout << "[WebP] requested_width=" << requested_width
                  << ", effective_width=" << BASE_IMAGE_WIDTH
                  << ", height=" << MAX_BLOCK_HEIGHT
                  << ", read_length=" << read_length
                  << ", sequences_per_row=" << sequences_per_row
                  << ", max_seqs_per_block=" << max_seqs_per_block
                  << std::endl;

        const WebpEncodeParams webp_params{
            webp_lossless,
            webp_quality,
            webp_method
        };

        ThreadPool thread_pool(static_cast<size_t>(thread_n));
        std::atomic<size_t> blocks_completed{0};
        std::atomic<size_t> total_blocks_enqueued{0};
        std::mutex           cout_mutex;
        std::exception_ptr   global_exception = nullptr;
        std::mutex           exception_mutex;

        // ★ 信号量：控制同时在飞的文件数
        std::atomic<int> inflight_files{0};
        std::mutex       inflight_mutex;
        std::condition_variable inflight_cv;

        for (int id = 0; id < TOTAL_FILES; ++id)
        {
            // ★ 等待 inflight 文件数低于上限
            {
                std::unique_lock<std::mutex> lk(inflight_mutex);
                inflight_cv.wait(lk, [&] {
                    return inflight_files.load() < MAX_INFLIGHT_FILES;
                });
                inflight_files++;
            }

            // 检查是否有异常
            {
                std::lock_guard<std::mutex> lk(exception_mutex);
                if (global_exception)
                    std::rethrow_exception(global_exception);
            }

            std::vector<std::string> raw_blocks;
            try {
                raw_blocks = ReadQualityScoreBlocks_40(id, quality_score_dir, read_length);
            } catch (const std::exception& e) {
                inflight_files--;
                inflight_cv.notify_one();
                std::lock_guard<std::mutex> lk(cout_mutex);
                std::cerr << "[警告] 读取文件" << id << "失败：" << e.what()
                          << "，跳过" << std::endl;
                continue;
            }

            if (raw_blocks.empty()) {
                inflight_files--;
                inflight_cv.notify_one();
                std::lock_guard<std::mutex> lk(cout_mutex);
                std::cout << "[警告] 文件" << id << "无有效质量块，跳过" << std::endl;
                continue;
            }

            auto matrix = std::make_shared<std::vector<std::vector<uint8_t>>>();
            matrix->reserve(raw_blocks.size());
            for (const auto& b : raw_blocks)
                matrix->emplace_back(b.begin(), b.end());

            const std::string output_path =
                std::string(webp_dir) + "/combined_output." + std::to_string(id);
            std::error_code mkdir_error;
            fs::create_directories(output_path, mkdir_error);
            if (mkdir_error) {
                throw std::runtime_error("创建目录失败: " + output_path +
                                         ": " + mkdir_error.message());
            }

            const size_t total_seqs = matrix->size();
            {
                std::ofstream count_file(output_path + "/sequence_count.txt",
                                         std::ios::out | std::ios::trunc);
                if (!count_file || !(count_file << total_seqs << '\n')) {
                    throw std::runtime_error("写入序列数量失败: " + output_path);
                }
            }
            const size_t num_blocks =
                (total_seqs + max_seqs_per_block - 1) / max_seqs_per_block;

            // ★ 用 shared_ptr 包装一个计数器，所有 block 完成后自动释放 inflight 槽
            auto remaining = std::make_shared<std::atomic<size_t>>(num_blocks);

            total_blocks_enqueued += num_blocks;
            if (id % 100 == 0 || id + 1 == TOTAL_FILES) {
                std::lock_guard<std::mutex> lk(cout_mutex);
                std::cout << "[入队] 文件" << id
                          << "，序列数=" << total_seqs
                          << "，WebP块数=" << num_blocks << std::endl;
            }

            for (size_t blk = 0; blk < num_blocks; ++blk) {
                const size_t start_seq = blk * max_seqs_per_block;
                const size_t end_seq   = std::min(start_seq + max_seqs_per_block, total_seqs);
                thread_pool.enqueue(
                    [matrix, output_path, read_length,
                     blk, start_seq, end_seq,
                     BASE_IMAGE_WIDTH, sequences_per_row,
                     webp_params,
                     &blocks_completed,
                     &global_exception, &exception_mutex,
                     remaining, &inflight_files, &inflight_cv]()
                    {
                        try {
                            EncodeSingleWebpBlock(
                                *matrix, output_path, read_length,
                                blk, start_seq, end_seq,
                                BASE_IMAGE_WIDTH, sequences_per_row,
                                webp_params);
                            ++blocks_completed;
                        } catch (...) {
                            std::lock_guard<std::mutex> lk(exception_mutex);
                            if (!global_exception)
                                global_exception = std::current_exception();
                            ++blocks_completed;
                        }

                        // ★ 这个文件最后一个 block 完成时释放 inflight 槽
                        if (--(*remaining) == 0) {
                            inflight_files--;
                            inflight_cv.notify_one();
                        }
                    });
            }
        }

        // 等待所有 block 完成
        const size_t expected = total_blocks_enqueued.load();
        const auto   start_time = std::chrono::steady_clock::now();
        while (blocks_completed.load() < expected) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            {
                std::lock_guard<std::mutex> lk(exception_mutex);
                if (global_exception)
                    std::rethrow_exception(global_exception);
            }
            const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start_time).count();
            std::lock_guard<std::mutex> lk(cout_mutex);
            std::cout << "全局进度: " << blocks_completed.load() << "/" << expected
                      << " 块完成 (耗时: " << elapsed << "秒，inflight文件="
                      << inflight_files.load() << ")" << std::endl;
        }

        {
            std::lock_guard<std::mutex> lk(exception_mutex);
            if (global_exception)
                std::rethrow_exception(global_exception);
        }
        std::cout << "\n所有块处理完成！共 " << expected
                  << " 块，read_length=" << read_length << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\n致命错误: " << e.what() << std::endl;
    }
}
