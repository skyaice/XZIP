#include <vector>
#include <string>
#include <webp/decode.h>
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <cstdint>
#include <dirent.h>
#include <sys/stat.h>
#include <algorithm>
#include <cstring>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <condition_variable>
#include <chrono>
#include <iomanip>
#include <functional>
#include <ctime>

const std::unordered_map<char, std::string> kCharToBinary_40 = {
    {'!', "000000"}, {'\"', "000001"}, {'#', "000010"}, {'$', "000011"},
    {'%', "000100"}, {'&', "000101"}, {'\'', "000110"}, {'(', "000111"},
    {')', "001000"}, {'*', "001001"}, {'+', "001010"}, {',', "001011"},
    {'-', "001100"}, {'.', "001101"}, {'/', "001110"}, {'0', "001111"},
    {'1', "010000"}, {'2', "010001"}, {'3', "010010"}, {'4', "010011"},
    {'5', "010100"}, {'6', "010101"}, {'7', "010110"}, {'8', "010111"},
    {'9', "011000"}, {':', "011001"}, {';', "011010"}, {'<', "011011"},
    {'=', "011100"}, {'>', "011101"}, {'?', "011110"}, {'@', "011111"},
    {'A', "100000"}, {'B', "100001"}, {'C', "100010"}, {'D', "100011"},
    {'E', "100100"}, {'F', "100101"}, {'G', "100110"}, {'H', "100111"}
};

const std::unordered_map<std::string, char> kBinaryToChar_40 = {
    {"000000", '!'}, {"000001", '\"'}, {"000010", '#'}, {"000011", '$'},
    {"000100", '%'}, {"000101", '&'}, {"000110", '\''}, {"000111", '('},
    {"001000", ')'}, {"001001", '*'}, {"001010", '+'}, {"001011", ','},
    {"001100", '-'}, {"001101", '.'}, {"001110", '/'}, {"001111", '0'},
    {"010000", '1'}, {"010001", '2'}, {"010010", '3'}, {"010011", '4'},
    {"010100", '5'}, {"010101", '6'}, {"010110", '7'}, {"010111", '8'},
    {"011000", '9'}, {"011001", ':'}, {"011010", ';'}, {"011011", '<'},
    {"011100", '='}, {"011101", '>'}, {"011110", '?'}, {"011111", '@'},
    {"100000", 'A'}, {"100001", 'B'}, {"100010", 'C'}, {"100011", 'D'},
    {"100100", 'E'}, {"100101", 'F'}, {"100110", 'G'}, {"100111", 'H'}
};

// 进度跟踪器结构体
struct ProgressTracker_40 {
    std::atomic<int> processed{0};
    std::atomic<int> total{0};
    std::atomic<int> success{0};
    std::atomic<int> errors{0};
    std::chrono::time_point<std::chrono::steady_clock> start_time;
};

// 创建多级目录
void CreateDirectoryRecursive_40(const char* path) {
    char tmp[1024];
    char* p = NULL;
    snprintf(tmp, sizeof(tmp), "%s", path);
    size_t len = strlen(tmp);
    if (tmp[len - 1] == '/') tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;
            mkdir(tmp, S_IRWXU);
            *p = '/';
        }
    }
    mkdir(tmp, S_IRWXU);
}

// 获取目录下的webp文件列表
void GetWebpFiles_40(const std::string& dir_path, std::vector<std::string>& files) {
    DIR* dir = opendir(dir_path.c_str());
    if (!dir) throw std::runtime_error("无法打开目录: " + dir_path);
    struct dirent* entry;
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_type == DT_REG) {
            const std::string name(entry->d_name);
            if (name.size() > 5 && name.compare(name.size() - 5, 5, ".webp") == 0) {
                files.push_back(dir_path + "/" + name);
            }
        }
    }
    closedir(dir);
}

// 自定义排序函数
bool WebpFileComparator_40(const std::string& a, const std::string& b) {
    size_t pos_a = a.find_last_of("_");
    size_t dot_a = a.find_last_of(".");
    size_t pos_b = b.find_last_of("_");
    size_t dot_b = b.find_last_of(".");
    int num_a = std::stoi(a.substr(pos_a + 1, dot_a - pos_a - 1));
    int num_b = std::stoi(b.substr(pos_b + 1, dot_b - pos_b - 1));
    return num_a < num_b;
}

std::vector<std::vector<uint8_t>> RestoreMatrixFromWebp_40(const std::string& input_dir, int read_length) {
    std::vector<std::vector<uint8_t>> restored_matrix;

    // 检查目录是否存在
    struct stat info;
    if (stat(input_dir.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        std::cerr << "警告: 目录不存在，跳过处理: " << input_dir << std::endl;
        return restored_matrix;
    }

    // 获取并排序文件列表
    std::vector<std::string> webp_files;
    try {
        GetWebpFiles_40(input_dir, webp_files);
        if (webp_files.empty()) {
            std::cerr << "警告: 目录中未找到WebP文件: " << input_dir << std::endl;
            return restored_matrix;
        }
        std::sort(webp_files.begin(), webp_files.end(), WebpFileComparator_40);
    } catch (const std::runtime_error& e) {
        std::cerr << "警告: " << e.what() << "，跳过目录: " << input_dir << std::endl;
        return restored_matrix;
    }

    // 严格对齐编码端的核心参数计算逻辑
    const int pixels_per_sequence = read_length;
    const int max_webp_dimension  = 4096;
    const int sequences_per_row   = max_webp_dimension / pixels_per_sequence;
    if (sequences_per_row == 0) {
        throw std::runtime_error("read_length=" + std::to_string(read_length) +
                                 " 超过WebP最大宽度限制(16383)");
    }
    const int BASE_IMAGE_WIDTH = sequences_per_row * pixels_per_sequence;

    // 处理每个WebP文件（按块顺序）
    for (const auto& file_path : webp_files) {
        // 1. 读取WebP文件数据
        std::ifstream file(file_path, std::ios::binary | std::ios::ate);
        if (!file.is_open()) {
            throw std::runtime_error("无法打开文件: " + file_path + " (" + strerror(errno) + ")");
        }
        size_t file_size = file.tellg();
        if (file_size == 0) {
            file.close();
            throw std::runtime_error("文件为空: " + file_path);
        }
        file.seekg(0, std::ios::beg);
        std::vector<uint8_t> webp_data(file_size);
        if (!file.read(reinterpret_cast<char*>(webp_data.data()), file_size)) {
            file.close();
            throw std::runtime_error("读取文件失败: " + file_path + " (" + strerror(errno) + ")");
        }
        file.close();

        // 2. 获取文件头信息（宽高）
        WebPDecoderConfig config;
        if (!WebPInitDecoderConfig(&config)) {
            throw std::runtime_error("WebP解码器初始化失败");
        }
        VP8StatusCode status = WebPGetFeatures(webp_data.data(), webp_data.size(), &config.input);
        if (status != VP8_STATUS_OK) {
            WebPFreeDecBuffer(&config.output);
            throw std::runtime_error("解析WebP头信息失败: " + file_path +
                                     " (错误码: " + std::to_string(status) + ")");
        }

        const int actual_width  = config.input.width;
        const int actual_height = config.input.height;

        // 允许 BASE_IMAGE_WIDTH 或 BASE_IMAGE_WIDTH+1（兼容偶数对齐变体）
        if (actual_width != BASE_IMAGE_WIDTH && actual_width != BASE_IMAGE_WIDTH + 1) {
            WebPFreeDecBuffer(&config.output);
            throw std::runtime_error(
                "WebP文件宽度不匹配: " + file_path +
                " 预期: " + std::to_string(BASE_IMAGE_WIDTH) +
                " 或 " + std::to_string(BASE_IMAGE_WIDTH + 1) +
                " 实际: " + std::to_string(actual_width));
        }

        // 3. RGB模式解码（无损ARGB编码，R通道值严格保留，alpha=255无预乘问题）
        int decoded_width = 0, decoded_height = 0;
        std::cout << "准备解码(RGB模式): " << file_path << "\n";
        uint8_t* rgb_data = WebPDecodeRGB(webp_data.data(), webp_data.size(),
                                          &decoded_width, &decoded_height);
        if (rgb_data == nullptr) {
            WebPFreeDecBuffer(&config.output);
            throw std::runtime_error("RGB模式解码失败: " + file_path);
        }
        std::cout << "成功解码(RGB模式): " << decoded_width << "x" << decoded_height << "\n";

        if (decoded_width != actual_width || decoded_height != actual_height) {
            WebPFree(rgb_data);
            WebPFreeDecBuffer(&config.output);
            throw std::runtime_error(
                "解码尺寸不匹配: " + file_path +
                " 头信息: " + std::to_string(actual_width) + "x" + std::to_string(actual_height) +
                " 解码结果: " + std::to_string(decoded_width) + "x" + std::to_string(decoded_height));
        }

        // 4. 提取矩阵数据
        //    【关键修复】编码端对图像块最后一行未填满的位置初始化为0，
        //    解码时需过滤这些全零填充序列。
        //    有效质量分数字符ASCII范围为33('!')~72('H')，永远不会为0，过滤安全。
        const int effective_width = std::min(actual_width, BASE_IMAGE_WIDTH);
        const int stride_rgb      = decoded_width * 3;

        size_t valid_seq_count = 0;

        for (int y = 0; y < decoded_height; ++y) {
            for (int seq_idx_in_row = 0; seq_idx_in_row < sequences_per_row; ++seq_idx_in_row) {

                const int x_start = seq_idx_in_row * pixels_per_sequence;
                std::vector<uint8_t> sequence(pixels_per_sequence, 0);
                bool seq_valid = true;

                for (int px = 0; px < pixels_per_sequence; ++px) {
                    const int x = x_start + px;
                    if (x >= effective_width) {
                        // 超出有效宽度，对应编码端未使用的填充列
                        sequence[px] = 0;
                        continue;
                    }
                    if (x >= decoded_width || y >= decoded_height) {
                        seq_valid = false;
                        break;
                    }
                    const int pixel_offset = y * stride_rgb + x * 3;
                    const size_t total_rgb_size = static_cast<size_t>(stride_rgb) * decoded_height;
                    if (pixel_offset + 2 >= static_cast<int>(total_rgb_size)) {
                        seq_valid = false;
                        break;
                    }
                    // 取R通道（编码端R=G=B=ascii_val，取任意一个均等价）
                    sequence[px] = rgb_data[pixel_offset];
                }

                if (seq_valid) {
                    // 【核心修复】过滤编码端零填充序列：
                    // 编码端对图像最后一行空余位置填充全0，
                    // 而有效质量分数ASCII >= 33，不可能出现0，直接过滤首字节为0的序列。
                    if (sequence[0] != 0) {
                        restored_matrix.push_back(std::move(sequence));
                        valid_seq_count++;
                    }
                }
            }
        }

        // 5. 资源清理
        WebPFree(rgb_data);
        WebPFreeDecBuffer(&config.output);

        std::cout << "已处理WebP文件: " << file_path
                  << " (尺寸: " << actual_width << "x" << actual_height
                  << ", 还原有效序列数: " << valid_seq_count << ")" << std::endl;
    }

    // 6. 最终验证：补全不足长度的序列（理论上不应出现，作为保险）
    if (!restored_matrix.empty()) {
        for (auto& seq : restored_matrix) {
            if (seq.size() != static_cast<size_t>(read_length)) {
                seq.resize(read_length, 0);
            }
        }
        if (restored_matrix[0].size() != static_cast<size_t>(read_length)) {
            throw std::runtime_error(
                "还原矩阵列数错误: 预期 " + std::to_string(read_length) +
                " 实际 " + std::to_string(restored_matrix[0].size()));
        }
    }

    std::cout << "矩阵还原完成，总有效序列数: " << restored_matrix.size() << std::endl;
    return restored_matrix;
}

bool CompareMatrices_40(const std::vector<std::vector<uint8_t>>& original,
                        const std::vector<std::vector<uint8_t>>& restored) {
    if (!restored.size()) return true;
    const size_t CHECK_LIMIT = original.size();

    for (size_t i = 0; i < CHECK_LIMIT; ++i) {
        if (i >= original.size() || i >= restored.size()) {
            std::cerr << "矩阵行数不足" << CHECK_LIMIT
                      << ": original=" << original.size()
                      << " restored=" << restored.size() << std::endl;
            return false;
        }
        if (original[i].size() != restored[i].size()) {
            std::cerr << "第 " << i << " 行列数不一致: "
                      << original[i].size() << " vs " << restored[i].size() << std::endl;
            return false;
        }
        for (size_t j = 0; j < original[i].size(); ++j) {
            if (original[i][j] != restored[i][j]) {
                std::cerr << "差异位置 (" << i << "," << j << "): "
                          << static_cast<int>(original[i][j])
                          << " vs " << static_cast<int>(restored[i][j]) << std::endl;
                return false;
            }
        }
    }

    if (original.size() > CHECK_LIMIT || restored.size() > CHECK_LIMIT) {
        std::cout << "前" << CHECK_LIMIT << "行验证通过（总行数: original="
                  << original.size() << " restored=" << restored.size() << "）" << std::endl;
    }
    return true;
}

// 线程安全队列
template <typename T>
class ThreadSafeQueue_40 {
public:
    void Push_40(T value) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push(std::move(value));
        cond_.notify_one();
    }

    bool Pop_40(T& value) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_.wait(lock, [this] { return !queue_.empty() || stop_; });
        if (stop_ && queue_.empty()) return false;
        value = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void Stop_40() {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
        cond_.notify_all();
    }

    bool Empty_40() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }

private:
    mutable std::mutex mutex_;
    std::condition_variable cond_;
    std::queue<T> queue_;
    bool stop_ = false;
};

// 处理单个ID的任务
void ProcessOneID_40(int id, const std::string& webp_base_dir,
                     const std::string& output_base_dir,
                     std::mutex& output_mutex,
                     std::atomic<bool>& error_occurred,
                     ProgressTracker_40& progress,
                     int read_length) {
    try {
        if (error_occurred.load()) return;

        std::string webp_dir = webp_base_dir + "combined_output." + std::to_string(id);
        std::cout << "正在准备还原Webp\n";
        auto quality_blocks = RestoreMatrixFromWebp_40(webp_dir, read_length);
        std::cout << "已经还原Webp\n";

        if (quality_blocks.empty()) {
            progress.errors++;
            progress.processed++;
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "跳过缺失目录: ID " << id << std::endl;
            return;
        }

        std::string output_file = output_base_dir + "quality_score." + std::to_string(id) + ".bin";
        {
            std::ofstream out(output_file, std::ios::binary);
            size_t block_count = quality_blocks.size();
            out.write(reinterpret_cast<const char*>(&block_count), sizeof(size_t));
            std::cout << output_file << " size: " << block_count << std::endl;
            for (const auto& block : quality_blocks) {
                size_t len = block.size();
                out.write(reinterpret_cast<const char*>(&len), sizeof(size_t));
                out.write(reinterpret_cast<const char*>(block.data()), len);
            }
        }
    } catch (const std::exception& e) {
        error_occurred = true;
        progress.errors++;
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "\n处理ID " << id << " 时出错: " << e.what() << std::endl;
    }
    progress.processed++;
}

// 工作线程函数
void WorkerThread_40(ThreadSafeQueue_40<int>& task_queue,
                     const std::string& webp_base_dir,
                     const std::string& output_base_dir,
                     std::mutex& output_mutex,
                     std::atomic<bool>& error_occurred,
                     ProgressTracker_40& progress,
                     int read_length) {
    int id;
    while (task_queue.Pop_40(id)) {
        if (error_occurred.load()) {
            task_queue.Stop_40();
            break;
        }
        ProcessOneID_40(id, webp_base_dir, output_base_dir,
                        output_mutex, error_occurred, progress, read_length);
    }
}

// 更新进度显示
void UpdateProgressDisplay_40(const ProgressTracker_40& progress) {
    auto now = std::chrono::steady_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now - progress.start_time);
    auto elapsed_millis  = std::chrono::duration_cast<std::chrono::milliseconds>(now - progress.start_time);

    int processed = progress.processed.load();
    int total     = progress.total.load();
    int success   = progress.success.load();
    int errors    = progress.errors.load();

    if (total <= 0) return;

    float percentage  = 100.0f * processed / total;
    const int bar_width = 50;
    int pos = static_cast<int>(bar_width * percentage / 100.0);

    std::cout << "\r[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] ";
    std::cout << std::setw(5) << std::fixed << std::setprecision(1) << percentage << "% "
              << "处理: " << processed << "/" << total << " "
              << "成功: " << success << " "
              << "错误: " << errors << " "
              << "用时: " << elapsed_seconds.count() << "秒";

    if (processed > 0 && processed < total) {
        auto remaining = std::chrono::duration_cast<std::chrono::seconds>(
            (now - progress.start_time) * (total - processed) / processed);
        std::cout << " 剩余: " << remaining.count() << "秒";
    }
    std::cout << std::flush;
}

int webp_reconstructor_40_main(char* webp_dir, char* output_dir,
                                int thread_n, int read_length, int total) {
    try {
        const int START_ID   = 0;
        const int END_ID     = total;
        const int NUM_THREADS = thread_n;

        std::string webp_base_dir = webp_dir;
        if (!webp_base_dir.empty() && webp_base_dir.back() != '/') webp_base_dir += '/';

        std::string output_base_dir = output_dir;
        if (!output_base_dir.empty() && output_base_dir.back() != '/') output_base_dir += '/';

        CreateDirectoryRecursive_40(output_base_dir.c_str());

        ProgressTracker_40 progress;
        progress.total      = END_ID - START_ID + 1;
        progress.start_time = std::chrono::steady_clock::now();

        std::time_t start_time_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::cout << "使用 " << NUM_THREADS << " 个线程进行解压处理" << std::endl;
        std::cout << "处理ID范围: [" << START_ID << ", " << END_ID << "]" << std::endl;
        std::cout << "开始时间: " << std::ctime(&start_time_t);
        std::cout << "源目录: " << webp_base_dir << std::endl;
        std::cout << "输出目录: " << output_base_dir << std::endl;

        ThreadSafeQueue_40<int> task_queue;
        for (int id = START_ID; id <= END_ID; ++id) task_queue.Push_40(id);
        task_queue.Stop_40();

        std::mutex          output_mutex;
        std::atomic<bool>   error_occurred(false);
        std::vector<std::thread> threads;
        std::cout << "work here\n";

        for (int i = 0; i < NUM_THREADS; ++i) {
            threads.emplace_back([&] {
                WorkerThread_40(task_queue, webp_base_dir, output_base_dir,
                                output_mutex, error_occurred, progress, read_length);
            });
        }

        while (!task_queue.Empty_40() && !error_occurred.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            std::lock_guard<std::mutex> lock(output_mutex);
            UpdateProgressDisplay_40(progress);
        }

        for (auto& thread : threads)
            if (thread.joinable()) thread.join();

        UpdateProgressDisplay_40(progress);
        std::cout << std::endl;

        if (error_occurred.load())
            throw std::runtime_error("处理过程中出现错误");

        auto end_time      = std::chrono::steady_clock::now();
        auto total_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - progress.start_time);
        auto total_millis  = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - progress.start_time);
        std::time_t end_time_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

        std::cout << "\n==========================================" << std::endl;
        std::cout << "        所有文件处理完成"                      << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "总任务数: " << progress.total   << std::endl;
        std::cout << "成功处理: " << progress.success << std::endl;
        std::cout << "失败任务: " << progress.errors  << std::endl;
        std::cout << "开始时间: " << std::ctime(&start_time_t);
        std::cout << "结束时间: " << std::ctime(&end_time_t);
        std::cout << "总耗时  : " << total_seconds.count() << " 秒 ("
                  << total_millis.count() << " 毫秒)" << std::endl;
        std::cout << "处理速度: " << std::fixed << std::setprecision(2)
                  << (progress.total * 1000.0 / total_millis.count())
                  << " 任务/秒" << std::endl;
        std::cout << "==========================================" << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n程序错误: " << e.what() << std::endl;
        return 1;
    }
}