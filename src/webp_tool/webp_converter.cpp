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
#include <sys/stat.h>
#include "../../../src/BWT_aln.hpp"

namespace fs = std::filesystem;

const std::unordered_map<char, std::string> kCharToBinary = {
    {'F', "00"},
    {':', "01"},
    {'#', "10"},
    {',', "11"}
};

constexpr uint8_t kBitMask[8] = {128, 64, 32, 16, 8, 4, 2, 1};

std::string ConvertToBinaryString(const std::string& quality_str, int read_length) {
    std::string binary;
    binary.reserve(2*read_length);
    
    for (char c : quality_str) {
        auto it = kCharToBinary.find(c);
        if (it == kCharToBinary.end()) {
            throw std::runtime_error("Invalid quality character: " + std::string(1, c));
        }
        binary += it->second;
    }
    
    if (binary.length() != 2*read_length) {
        throw std::runtime_error("Binary conversion error");
    }
    return binary;
}

std::vector<uint8_t> BinaryToBytes(const std::string& binary, int read_length) {
    // 计算目标总位数和字节数：2*read_length位 → 向上取整为字节数
    const int total_bits = 2 * read_length;
    const int byte_count = (total_bits + 7) / 8; // 等价于 ceil(total_bits/8)，避免截断
    
    // 预分配内存（无扩容）
    std::vector<uint8_t> bytes(byte_count, 0);
    uint8_t* byte_ptr = bytes.data();

    // 遍历二进制字符串，每8位转换为1个字节
    for (int i = 0; i < total_bits; i += 8) {
        uint8_t current_byte = 0;
        // 处理当前字节的8位（不足8位则补0）
        for (int j = 0; j < 8; ++j) {
            const int bit_pos = i + j;
            // 仅当位位置有效且为'1'时，设置对应位
            if (bit_pos < static_cast<int>(binary.size()) && binary[bit_pos] == '1') {
                current_byte |= kBitMask[j]; // 位或运算设置对应位
            }
            // 若bit_pos超出binary长度，默认补0（bytes初始化已为0，无需处理）
        }
        *byte_ptr++ = current_byte;
    }

    return bytes;
}


class ThreadPool {
    public:
        ThreadPool(size_t threads) : stop(false) {
            for(size_t i = 0; i < threads; ++i)
                workers.emplace_back([this] {
                    for(;;) {
                        std::function<void()> task;
                        {
                            std::unique_lock<std::mutex> lock(this->queue_mutex);
                            this->condition.wait(lock,
                                [this]{ return this->stop || !this->tasks.empty(); });
                            if(this->stop && this->tasks.empty())
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
            for(std::thread &worker: workers)
                worker.join();
        }
    
    private:
        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;
        std::mutex queue_mutex;
        std::condition_variable condition;
        bool stop;
};

void quality_diff_base_build(std::vector<std::string>quality_diff_base_blocks, const std::string quality_diff_base)
{
    std::unordered_map<char, uint64_t> freqMap;
    for (std::string& seq : quality_diff_base_blocks) {
        for (char c : seq) {
            freqMap[c]++;
        }
        freqMap['\0']++; // 为每个序列添加结束符
    }
    if (freqMap.empty()) {
        std::cerr << "diff_seq_build: 没有统计到任何字符，无法构建哈夫曼树。" << std::endl;
        return;
    }
    BWT_aln::HuffmanCoder coder;
    auto root = coder.buildHuffmanTree(freqMap);
    std::unordered_map<char, std::string> huffmanCode;
    coder.buildCodes(root, "", huffmanCode);
    std::string full_binary_string;
    for (const std::string& seq : quality_diff_base_blocks) {
        for (char c : seq) {
            full_binary_string += huffmanCode[c];
        }
        full_binary_string += huffmanCode['\0']; // 添加序列结束符的编码
    }
    if (full_binary_string.empty()) {
        std::cerr << "quality_diff_base_build: 编码后没有生成任何二进制数据，无法写入文件。" << std::endl;
        return;
    }
    
    // 保存原始二进制字符串的长度，用于解码时截断
    uint64_t original_bits_length = full_binary_string.size();

    // 4. 将01字符串打包成uint8_t字节数组
    std::vector<uint8_t> byte_data;
    for (size_t i = 0; i < full_binary_string.size(); i += 8) {
        std::string byte_str = full_binary_string.substr(i, 8);
        while (byte_str.size() < 8) {
            byte_str += '0'; // 不足8位补0
        }
        uint8_t byte = static_cast<uint8_t>(std::stoi(byte_str, nullptr, 2));
        byte_data.push_back(byte);
    }

    std::ofstream ofs(quality_diff_base, std::ios::binary);
    if (!ofs.is_open()) {
        std::cerr << "diff_seq_build: 无法打开文件 " << quality_diff_base << " 进行写入。" << std::endl;
        return;
    }

    const char magic[4] = { 'H', 'u', 'f', 'f' };
    ofs.write(magic, sizeof(magic));

    uint8_t code_table_size = static_cast<uint8_t>(huffmanCode.size());
    ofs.write(reinterpret_cast<const char*>(&code_table_size), sizeof(code_table_size));

    // 3. 写入编码表条目
    for (const auto& pair : huffmanCode) {
        char c = pair.first;
        const std::string& code = pair.second;
        ofs.write(&c, sizeof(c));
        uint8_t code_len = static_cast<uint8_t>(code.size());
        ofs.write(reinterpret_cast<const char*>(&code_len), sizeof(code_len));
        ofs.write(code.c_str(), code_len);
    }
    ofs.write(reinterpret_cast<const char*>(&original_bits_length), sizeof(original_bits_length));
    ofs.write(reinterpret_cast<const char*>(byte_data.data()), byte_data.size() * sizeof(uint8_t));

    if (ofs.good()) {
        std::cout << "quality_diff_base_build: 成功将 " << (ofs.tellp()) << " 字节的数据写入到 " << quality_diff_base << std::endl;
    } else {
        std::cerr << "quality_diff_base_build: 写入文件 " << quality_diff_base << " 时发生错误。" << std::endl;
    }

    ofs.close();
}

void SaveMatrixToWebpLossless(const std::vector<std::vector<uint8_t>>& matrix,
                             const std::string& output_dir, int read_length) 
{
    if(!matrix.size())return;

    const int pixels_per_sequence = read_length / 4 + 1;
    // 参数校验（与原代码相同）
    if (matrix.empty() || matrix[0].size() != pixels_per_sequence) {
        throw std::invalid_argument("输入矩阵列数错误");
    }

    // 创建输出目录（与原代码相同）
    if (system(("mkdir -p " + output_dir).c_str()) != 0) {
        throw std::runtime_error("创建目录失败: " + output_dir);
    }

    // 核心参数（与原代码相同）
    const int sequences_per_row = 16384 / pixels_per_sequence;
    const int width = sequences_per_row * pixels_per_sequence;
    const int max_webp_height = width;
    const int stride = (width * 3 + 3) & ~3;
    const int total_sequences = matrix.size();

    // 无损压缩专用配置
    WebPConfig config;
    if (!WebPConfigInit(&config)) {
        throw std::runtime_error("WebP配置初始化失败");
    }
    config.lossless = 1;          // 启用无损模式
    config.method = 6;            // 最高压缩级别 (0-6)
    config.quality = 20;          // 在无损模式下控制压缩速度/质量

    // 分块处理控制（与原代码相同）
    int current_sequence = 0;
    int webp_file_index = 0;

    while (current_sequence < total_sequences) 
    {
        // 分块参数计算（与原代码相同）
        const int sequences_remaining = total_sequences - current_sequence;
        const int max_sequences_in_webp = sequences_per_row * max_webp_height;
        const int sequences_in_this_webp = std::min(max_sequences_in_webp, sequences_remaining);
        const int actual_height = (sequences_in_this_webp + sequences_per_row - 1) / sequences_per_row;

        // 准备图像数据（与原代码相同）
        std::vector<uint8_t> rgb_data(stride * actual_height, 0);
        for (int seq_offset = 0; seq_offset < sequences_in_this_webp; ++seq_offset) 
        {
            const int global_seq_idx = current_sequence + seq_offset;
            const int row_in_webp = seq_offset / sequences_per_row;
            const int col_in_row = seq_offset % sequences_per_row;
            const int x_start = col_in_row * pixels_per_sequence;
            
            for (int pixel = 0; pixel < pixels_per_sequence; ++pixel) {
                const int x = x_start + pixel;
                const int base_idx = row_in_webp * stride + x * 3;
                rgb_data[base_idx] = matrix[global_seq_idx][pixel];
            }
        }

        // 无损编码核心逻辑
        WebPPicture pic;
        if (!WebPPictureInit(&pic)) {
            throw std::runtime_error("WebP图片初始化失败");
        }
        pic.width = width;
        pic.height = actual_height;
        pic.use_argb = 1;  // 必须启用ARGB格式
        
        // 导入RGB数据
        if (!WebPPictureImportRGB(&pic, rgb_data.data(), stride)) {
            WebPPictureFree(&pic);
            throw std::runtime_error("数据导入失败");
        }

        // 内存写入器配置
        WebPMemoryWriter writer;
        WebPMemoryWriterInit(&writer);
        pic.writer = WebPMemoryWrite;
        pic.custom_ptr = &writer;

        // 执行编码
        if (!WebPEncode(&config, &pic)) {
            WebPMemoryWriterClear(&writer);
            WebPPictureFree(&pic);
            throw std::runtime_error("无损编码失败");
        }

        // 写入文件
        const std::string filename = output_dir + "/lossless_block_" + std::to_string(webp_file_index++) + ".webp";
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs.write(reinterpret_cast<char*>(writer.mem), writer.size)) {
            WebPMemoryWriterClear(&writer);
            WebPPictureFree(&pic);
            throw std::runtime_error("文件写入失败: " + filename);
        }

        // 清理资源
        WebPMemoryWriterClear(&writer);
        WebPPictureFree(&pic);

        // 更新进度（与原代码相同）
        current_sequence += sequences_in_this_webp;
        std::cout << "已生成无损文件: " << filename 
                  << " (序列 " << current_sequence - sequences_in_this_webp << "-"
                  << current_sequence - 1 << ")" << std::endl;
    }
}

void SaveAsciiMatrixToWebpLossless(const std::vector<std::vector<uint8_t>>& matrix,
                                   const std::string& output_dir_path, int read_length) // 现在接收的是文件夹路径
{
    if (matrix.empty()) {
        std::cout << "输入矩阵为空，跳过WebP生成。" << std::endl;
        return;
    }

    // 1. 核心配置参数
    const int pixels_per_sequence = read_length;  // 每个质量串占read_length个像素
    const int sequences_per_row = 16384 / pixels_per_sequence;    // 每行放置质量串
    const int image_width = sequences_per_row * pixels_per_sequence; // 图像宽度 = 109 * 150
    const int max_webp_dimension = 16383; // WebP最大单维度尺寸 (修正为16383)

    // 2. 验证输入矩阵的列数是否正确
    if (matrix[0].size() != pixels_per_sequence) {
        throw std::invalid_argument("输入矩阵的列数必须为 " + std::to_string(pixels_per_sequence) +
                                    "，当前为 " + std::to_string(matrix[0].size()));
    }

    // 3. WebP编码配置
    WebPConfig config;
    if (!WebPConfigInit(&config)) {
        throw std::runtime_error("WebP配置初始化失败");
    }
    config.lossless = 1; // 启用无损压缩
    config.method = 6;   // 最高压缩级别 (0-6)

    // 4. 计算分块信息
    const int total_sequences = static_cast<int>(matrix.size());
    // 计算每个块最多能容纳的序列数
    const int max_sequences_per_block = sequences_per_row * max_webp_dimension;
    const int num_blocks = (total_sequences + max_sequences_per_block - 1) / max_sequences_per_block;

    std::cout << "开始生成WebP文件，总共 " << total_sequences << " 个序列，将分为 " << num_blocks << " 个文件。" << std::endl;

    // 5. 循环处理每个块
    for (int block_idx = 0; block_idx < num_blocks; ++block_idx) {
        const int start_seq = block_idx * max_sequences_per_block;
        const int end_seq = std::min(start_seq + max_sequences_per_block, total_sequences);
        const int seqs_in_this_block = end_seq - start_seq;

        // 计算当前块的图像高度
        const int block_height = (seqs_in_this_block + sequences_per_row - 1) / sequences_per_row;

        // 检查尺寸是否超限
        if (block_height > max_webp_dimension) {
             std::cerr << "警告: 块 " << block_idx << " 的高度 " << block_height << " 超过了WebP限制 " << max_webp_dimension << "，可能导致编码失败。" << std::endl;
        }

        // 分配RGB数据缓冲区 (WebP需要RGB格式)
        const int stride = image_width * 3; // 每行的字节数 (RGB)
        std::vector<uint8_t> rgb_data(stride * block_height, 0);

        // 填充图像数据
        for (int seq_offset = 0; seq_offset < seqs_in_this_block; ++seq_offset) {
            const int global_seq_idx = start_seq + seq_offset;
            const int row_in_image = seq_offset / sequences_per_row;
            const int col_in_row = seq_offset % sequences_per_row;
            const int x_start = col_in_row * pixels_per_sequence;

            // 将ASCII码填充到RGB三个通道，形成灰度图像
            for (int pixel = 0; pixel < pixels_per_sequence; ++pixel) {
                const int x = x_start + pixel;
                const int base_idx = row_in_image * stride + x * 3;
                uint8_t ascii_value = matrix[global_seq_idx][pixel];
                rgb_data[base_idx] = ascii_value;     // R
                rgb_data[base_idx + 1] = ascii_value; // G
                rgb_data[base_idx + 2] = ascii_value; // B
            }
        }

        // 准备WebP图片结构
        WebPPicture pic;
        if (!WebPPictureInit(&pic)) {
            throw std::runtime_error("WebP图片初始化失败");
        }
        pic.width = image_width;
        pic.height = block_height;
        pic.use_argb = 0; // 使用RGB模式

        // 导入RGB数据
        if (!WebPPictureImportRGB(&pic, rgb_data.data(), stride)) {
            WebPPictureFree(&pic);
            throw std::runtime_error("WebP数据导入失败");
        }

        // 内存写入器设置
        WebPMemoryWriter writer;
        WebPMemoryWriterInit(&writer);
        pic.writer = WebPMemoryWrite;
        pic.custom_ptr = &writer;

        // 执行编码
        if (!WebPEncode(&config, &pic)) {
            WebPMemoryWriterClear(&writer);
            WebPPictureFree(&pic);
            throw std::runtime_error("WebP无损编码失败");
        }

        // --- 【核心修改点】构建输出文件路径 ---
        std::string filename = "part_" + std::to_string(block_idx) + ".webp";
        std::string file_path = output_dir_path + "/" + filename; // 在指定文件夹内创建文件

        // 写入文件
        std::ofstream ofs(file_path, std::ios::binary);
        if (!ofs.is_open()) {
            WebPMemoryWriterClear(&writer);
            WebPPictureFree(&pic);
            throw std::runtime_error("无法打开文件进行写入: " + file_path);
        }
        if (!ofs.write(reinterpret_cast<char*>(writer.mem), writer.size)) {
            WebPMemoryWriterClear(&writer);
            WebPPictureFree(&pic);
            throw std::runtime_error("文件写入失败: " + file_path);
        }
        ofs.close();

        // 清理资源
        WebPMemoryWriterClear(&writer);
        WebPPictureFree(&pic);

        std::cout << "已生成WebP文件: " << file_path
                  << " (包含序列 " << start_seq << " 到 " << end_seq - 1 << ")" << std::endl;
    }
}


void ProcessQualityBlocks(const std::vector<std::string>& quality_blocks, 
                         const std::string& output_prefix, int file_id, int read_length, const std::string quality_diff_base) {
    std::vector<std::vector<uint8_t>> image_matrix;
    /*
    for (const auto& block : quality_blocks) {
        if (block.length() != read_length) {
            throw std::runtime_error("Invalid block length: " + std::to_string(block.length()));
        }
        
        std::string binary = ConvertToBinaryString(block, read_length);
        image_matrix.push_back(BinaryToBytes(binary, read_length));
    }
    */
    std::vector<std::string>quality_diff_base_blocks;
    for (const auto& block : quality_blocks) {
        if (block.length() != 150) {
            throw std::runtime_error("无效的块长度: " + std::to_string(block.length()) +
                                     "。每个块必须是150个字符。");
        }
        std::vector<uint8_t> ascii_row;
        std::string diff_base;
        ascii_row.reserve(150);
        diff_base.clear();
        for (char c : block) {
            if(c=='F')ascii_row.push_back(0);
            else
            {
                ascii_row.push_back(255);
                diff_base += c;
            }
        }
        image_matrix.push_back(ascii_row);
        quality_diff_base_blocks.push_back(diff_base);
    }
    if(image_matrix.size())std::cout<<"matrix size: "<<image_matrix.size()<<" "<<image_matrix[0].size()<<std::endl;
    //SaveMatrixToWebpLossless(image_matrix, output_prefix, read_length);
    SaveAsciiMatrixToWebpLossless(image_matrix, output_prefix, read_length);
    quality_diff_base_build(quality_diff_base_blocks, quality_diff_base);
}
std::vector<std::string> ReadQualityScoreBlocks(int id, char *quality_score_dir, int read_length) {
    std::string filename = std::string(quality_score_dir) + "/quality_score." + std::to_string(id) + ".bin";

    std::ifstream in_file(filename, std::ios::binary);
    if (!in_file) {
        throw std::runtime_error("cannot open the file: " + filename);
    }

    size_t block_count = 0;
    in_file.read(reinterpret_cast<char*>(&block_count), sizeof(size_t));

    std::vector<std::string> blocks;
    blocks.reserve(block_count);

    for (size_t i = 0; i < block_count; ++i) {
        size_t code_len = 0;
        in_file.read(reinterpret_cast<char*>(&code_len), sizeof(size_t));

        std::string str(code_len, '\0');
        in_file.read(&str[0], code_len);

        if (code_len != read_length) {
            std::cerr << "warning: file" << id << "the" << i 
                      << "block length is error: " << code_len << std::endl;
        }
        else blocks.emplace_back(std::move(str));
    }


    if (in_file.fail() && !in_file.eof()) {
        throw std::runtime_error("cannot read the file: " + filename);
    }

    return blocks;
}

void webp_main(char *quality_score_dir, char *webp_dir, int thread_n, int total, int read_length) {
    try {
        const int TOTAL_FILES = total + 1;
        const int START_ID = 0;
        const int NUM_FILES = TOTAL_FILES - START_ID;

        // 线程池配置
        const unsigned int num_workers = thread_n;
        ThreadPool thread_pool(num_workers);
        
        // 并发控制
        std::atomic<int> files_completed{0};
        std::mutex cout_mutex;
        std::exception_ptr global_exception = nullptr;
        std::mutex exception_mutex;

        // 提交任务
        for (int id = START_ID; id < TOTAL_FILES; ++id) {
            thread_pool.enqueue([id, webp_dir, quality_score_dir, &files_completed, &cout_mutex, 
                &global_exception, &exception_mutex, read_length]() 
            {
                try {
                    const auto blocks = ReadQualityScoreBlocks(id, quality_score_dir, read_length);
                    if(blocks.size()==0)
                    {
                        std::lock_guard<std::mutex> lock(cout_mutex);
                        std::cout << "[Warning] No blocks found in file " << id << std::endl;
                        ++files_completed;
                        return;
                    }
                    
                    // 生成输出路径
                    const std::string output_path = std::string(webp_dir) + "/combined_output." + std::to_string(id);
                    std::string dir(webp_dir);
                    size_t last_slash = dir.find_last_of("/\\");
                    std::string parent_dir;
                    if (last_slash == std::string::npos) {
                        parent_dir =  ".";
                    }
                    else parent_dir = dir.substr(0, last_slash);
                    const std::string quality_diff_base_dir = parent_dir + "/quality_diff_base";
                    mkdir(quality_diff_base_dir.c_str(), 0755);
                    const std::string quality_diff_base = quality_diff_base_dir + "/diff_" + std::to_string(id) + ".bin";
                    
                    // 创建输出目录（线程安全）
                    {
                        if (system(("mkdir -p " + output_path).c_str()) != 0) {
                            throw std::runtime_error("Failed to create directory: " + output_path);
                        }
                    }

                    // 处理并保存
                    ProcessQualityBlocks(blocks, output_path, id, read_length, quality_diff_base);

                    // 输出进度（线程安全）
                    {
                        std::cout << "[Success] Processed file " << id 
                                  << " (" << blocks.size() << " blocks)" 
                                  << std::endl;
                    }

                    ++files_completed;
                } catch (...) {
                    std::lock_guard<std::mutex> lock(exception_mutex);
                    if (!global_exception) {
                        global_exception = std::current_exception();
                    }
                }
            });
        }

        // 等待任务完成
        const auto start_time = std::chrono::steady_clock::now();
        while (files_completed < NUM_FILES) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            
            // 异常检查
            std::lock_guard<std::mutex> lock(exception_mutex);
            if (global_exception) {
                std::rethrow_exception(global_exception);
            }

            // 输出进度
            const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start_time
            ).count();
            
            std::cout << "Progress: " << files_completed << "/" << NUM_FILES
                      << " (" << elapsed << "s elapsed)" 
                      << std::endl;
        }

        std::cout << "\nAll files processed successfully. Total: " << NUM_FILES 
                  << std::endl;
        return;
    } catch (const std::exception& e) {
        std::cerr << "\nFatal error: " << e.what() << std::endl;
        return;
    }
}
    
