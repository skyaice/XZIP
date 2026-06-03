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
#include "../../../src/BWT_aln.hpp"

const std::unordered_map<char, std::string> kCharToBinary = {
    {'F', "00"}, {':', "01"}, {'#', "10"}, {',', "11"}
};

const std::unordered_map<std::string, char> kBinaryToChar = {
    {"00", 'F'}, {"01", ':'}, {"10", '#'}, {"11", ','}
};

// 创建多级目录
void CreateDirectoryRecursive(const char* path) {
    char tmp[1024];
    char* p = NULL;
    snprintf(tmp, sizeof(tmp), "%s", path);
    size_t len = strlen(tmp);
    if (tmp[len - 1] == '/') {
        tmp[len - 1] = 0;
    }
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
void GetWebpFiles(const std::string& dir_path, std::vector<std::string>& files) {
    std::cout<<"dir_path: "<<dir_path<<std::endl;
    DIR* dir = opendir(dir_path.c_str());
    if (!dir) {
        throw std::runtime_error("无法打开目录: " + dir_path);
    }
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

std::string get_webp_parent_dir(const char* pos_dir) {
    std::string dir(pos_dir);
    // 查找最后一个路径分隔符（兼容Linux和Windows）
    size_t last_slash = dir.find_last_of("/\\");
    if (last_slash == std::string::npos) {
        // 没有找到分隔符，父路径为当前目录（或根据需求返回空）
        return ".";
    }
    // 返回分隔符之前的部分（不包含分隔符）
    return dir.substr(0, last_slash);
}

// 自定义排序函数
bool WebpFileComparator(const std::string& a, const std::string& b) {
    size_t pos_a = a.find_last_of("_");
    size_t dot_a = a.find_last_of(".");
    size_t pos_b = b.find_last_of("_");
    size_t dot_b = b.find_last_of(".");
    int num_a = std::stoi(a.substr(pos_a + 1, dot_a - pos_a - 1));
    int num_b = std::stoi(b.substr(pos_b + 1, dot_b - pos_b - 1));
    return num_a < num_b;
}

std::vector<std::vector<uint8_t>> RestoreMatrixFromWebpTiled(const std::string& file_path, int pixels_per_sequence) {
    std::vector<std::vector<uint8_t>> restored_matrix;

    if (pixels_per_sequence <= 0) {
        throw std::invalid_argument("pixels_per_sequence 必须是正数。");
    }

    // 1. 读取文件内容到内存
    std::ifstream file(file_path, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开文件: " + file_path);
    }
    size_t size = file.tellg();
    file.seekg(0, std::ios::beg);
    std::vector<uint8_t> webp_data(size);
    if (!file.read(reinterpret_cast<char*>(webp_data.data()), size)) {
        throw std::runtime_error("读取文件失败: " + file_path);
    }
    file.close();

    // 2. 使用WebP库进行解码
    WebPDecoderConfig config;
    if (!WebPInitDecoderConfig(&config)) {
        throw std::runtime_error("WebP解码器初始化失败");
    }

    if (WebPGetFeatures(webp_data.data(), webp_data.size(), &config.input) != VP8_STATUS_OK) {
        WebPFreeDecBuffer(&config.output);
        throw std::runtime_error("解析WebP头信息失败: " + file_path);
    }

    config.output.colorspace = MODE_RGBA;
    if (WebPDecode(webp_data.data(), webp_data.size(), &config) != VP8_STATUS_OK) {
        WebPFreeDecBuffer(&config.output);
        throw std::runtime_error("解码失败: " + file_path);
    }

    const int img_width = config.input.width;
    const int img_height = config.input.height;
    const uint8_t* rgba_data = config.output.u.RGBA.rgba;
    const int stride = config.output.u.RGBA.stride;

    // 3. 逆向编码时的平铺逻辑，恢复原始矩阵
    for (int y = 0; y < img_height; ++y) { // 遍历图像的每一行
        // 遍历当前行，并按 pixels_per_sequence 的宽度切分
        for (int x = 0; x < img_width; x += pixels_per_sequence) {
            std::vector<uint8_t> original_sequence;
            original_sequence.reserve(pixels_per_sequence);

            // 切分出一个原始序列
            for (int i = 0; i < pixels_per_sequence; ++i) {
                // 边界检查，防止越界（虽然理论上不应该发生）
                if (x + i >= img_width) {
                    break; 
                }

                size_t pixel_idx = static_cast<size_t>(y) * static_cast<size_t>(stride) + static_cast<size_t>(x + i) * 4;
                uint8_t gray_value = rgba_data[pixel_idx]; // 取R通道的值

                // 根据编码规则，将非0值视为255
                original_sequence.push_back(gray_value == 0 ? 0 : 255);
            }
            
            // 只有当切分出的序列长度正确时才添加
            if (original_sequence.size() == pixels_per_sequence) {
                restored_matrix.push_back(original_sequence);
            }
        }
    }

    // 4. 释放资源
    WebPFreeDecBuffer(&config.output);

    return restored_matrix;
}

std::vector<std::vector<uint8_t>> RestoreMatrixFromWebp(const std::string& input_dir) {
    std::vector<std::vector<uint8_t>> restored_matrix;

    // 获取并排序文件列表
    std::vector<std::string> webp_files;
    GetWebpFiles(input_dir, webp_files);
    std::sort(webp_files.begin(), webp_files.end(), WebpFileComparator);

    // 处理每个文件
    for (const auto& file_path : webp_files) {
        std::ifstream file(file_path, std::ios::binary | std::ios::ate);
        size_t size = file.tellg();
        file.seekg(0, std::ios::beg);
        std::vector<uint8_t> webp_data(size);
        if (!file.read(reinterpret_cast<char*>(webp_data.data()), size)) {
            throw std::runtime_error("读取文件失败: " + file_path);
        }

        WebPDecoderConfig config;
        if (!WebPInitDecoderConfig(&config)) {
            throw std::runtime_error("WebP解码器初始化失败");
        }

        if (WebPGetFeatures(webp_data.data(), webp_data.size(), &config.input) != VP8_STATUS_OK) {
            throw std::runtime_error("解析WebP头信息失败: " + file_path);
        }

        config.output.colorspace = MODE_RGBA;
        if (WebPDecode(webp_data.data(), webp_data.size(), &config) != VP8_STATUS_OK) {
            WebPFreeDecBuffer(&config.output);
            throw std::runtime_error("解码失败: " + file_path);
        }

        const int width = config.input.width;
        const int height = config.input.height;
        const uint8_t* rgba = config.output.u.RGBA.rgba;
        const int stride = config.output.u.RGBA.stride;

        // 提取数据到矩阵
        for (int y = 0; y < height; ++y) {
            for (int seq_x = 0; seq_x < 431; ++seq_x) {
                std::vector<uint8_t> sequence(38);
                const int x_start = seq_x * 38;
                for (int px = 0; px < 38; ++px) {
                    const int pixel_idx = y * stride + (x_start + px) * 4;
                    sequence[px] = rgba[pixel_idx];
                }
                restored_matrix.push_back(sequence);
            }
        }
        WebPFreeDecBuffer(&config.output);
    }

    if (!restored_matrix.empty() && restored_matrix[0].size() != 38) {
        throw std::runtime_error("还原矩阵列数错误");
    }

    return restored_matrix;
}

std::string BytesToQualityString(const std::vector<uint8_t>& bytes) {
    if (bytes.size() != 38) {
        throw std::invalid_argument("字节向量长度必须为38");
    }

    std::string binary_str;
    binary_str.reserve(304);
    for (uint8_t byte : bytes) {
        for (int i = 7; i >= 0; --i) {
            binary_str += (byte & (1 << i)) ? '1' : '0';
        }
    }
    binary_str.resize(300); // 移除填充位

    std::string quality_str;
    for (size_t i = 0; i < 300; i += 2) {
        std::string pair = binary_str.substr(i, 2);
        auto it = kBinaryToChar.find(pair);
        if (it == kBinaryToChar.end()) {
            throw std::runtime_error("无效的二进制对: " + pair);
        }
        quality_str += it->second;
    }

    return quality_str;
}

bool CompareMatrices(const std::vector<std::vector<uint8_t>>& original, const std::vector<std::vector<uint8_t>>& restored) {
    if(!restored.size()) return true;
    const size_t CHECK_LIMIT = original.size(); // 设置需要检查的行数阈值
    
    // 检查前CHECK_LIMIT行数据
    for (size_t i = 0; i < CHECK_LIMIT; ++i) {
        // 边界检查
        if (i >= original.size() || i >= restored.size()) {
            std::cerr << "矩阵行数不足" << CHECK_LIMIT << ": original=" << original.size() << " restored=" << restored.size() << std::endl;
            return false;
        }
        
        // 列数一致性检查
        if (original[i].size() != restored[i].size()) {
            std::cerr << "第 " << i << " 行列数不一致: " << original[i].size() << " vs " << restored[i].size() << std::endl;
            return false;
        }
        
        // 元素级比较
        for (size_t j = 0; j < original[i].size(); ++j) {
            if (original[i][j] != restored[i][j]) {
                std::cerr << "差异位置 (" << i << "," << j << "): " << static_cast<int>(original[i][j]) << " vs " << static_cast<int>(restored[i][j]) << std::endl;
                return false;
            }
        }
    }
    
    // 如果原始数据超过检查阈值，输出提示信息
    if (original.size() > CHECK_LIMIT || restored.size() > CHECK_LIMIT) {
        std::cout << "前" << CHECK_LIMIT << "行验证通过（总行数: original=" << original.size() << " restored=" << restored.size() << "）" << std::endl;
    }
    
    return true;
}

// 线程安全队列
template <typename T>
class ThreadSafeQueue {
public:
    void Push(T value) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push(std::move(value));
        cond_.notify_one();
    }

    bool Pop(T& value) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_.wait(lock, [this]{ return !queue_.empty() || stop_; });
        if (stop_) return false;
        value = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void Stop() {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
        cond_.notify_all();
    }

    bool Empty() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }

private:
    mutable std::mutex mutex_;
    std::condition_variable cond_;
    std::queue<T> queue_;
    bool stop_ = false;
};

std::vector<std::string>get_webp_to_diff_seq_blocks(const std::string& webp_base_dir, int read_length, int id)
{
    std::vector<std::string> diff_seq_blocks;
    std::string webp_dir = webp_base_dir + "/combined_output." + std::to_string(id);
    std::vector<std::string> webp_files;
    DIR* dir = opendir(webp_dir.c_str());
    if (!dir) {
        std::cerr << "webp_read: 无法打开目录: " << webp_dir << std::endl;
        return diff_seq_blocks;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_type == DT_REG) {
            std::string filename(entry->d_name);
            if (filename.size() >= 5 && filename.compare(filename.size() - 5, 5, ".webp") == 0) {
                webp_files.push_back(webp_dir + "/" + filename);
            }
        }
    }
    closedir(dir);
    
    if (webp_files.empty()) {
        std::cerr << "webp_read: 在目录 " << webp_dir << " 中未找到任何 .webp 文件。" << std::endl;
        return diff_seq_blocks;
    }
    std::sort(webp_files.begin(), webp_files.end(), [](const std::string& a, const std::string& b) {
        // 提取文件名中的数字部分进行比较
        size_t pos_a = a.find_last_of('_');
        size_t dot_a = a.find_last_of('.');
        int num_a = std::stoi(a.substr(pos_a + 1, dot_a - pos_a - 1));
            
        size_t pos_b = b.find_last_of('_');
        size_t dot_b = b.find_last_of('.');
        int num_b = std::stoi(b.substr(pos_b + 1, dot_b - pos_b - 1));
            
        return num_a < num_b;
    });
    std::vector<std::vector<uint8_t>> complete_matrix;
    try {
        for (const std::string& file_path : webp_files) {
            std::cout << "webp_read: 正在解码文件: " << file_path << std::endl;
            auto partial_matrix = RestoreMatrixFromWebpTiled(file_path, read_length);
            complete_matrix.insert(complete_matrix.end(), partial_matrix.begin(), partial_matrix.end());
        }
    } catch (const std::exception& e) {
        std::cerr << "webp_read: 解码WebP文件时发生错误: " << e.what() << std::endl;
        return diff_seq_blocks;
    }
    if (complete_matrix.empty()) {
        std::cerr << "webp_read: 解码后未得到任何像素数据。" << std::endl;
        return diff_seq_blocks;
    }

    // 将像素矩阵转换回字符串向量
    for (const auto& pixel_row : complete_matrix) {
        std::string seq;
        seq.reserve(pixel_row.size());
        for (uint8_t pixel : pixel_row) {
            seq += (pixel == 0) ? '0' : '1';
        }
        if (!seq.empty()) {
            diff_seq_blocks.push_back(seq);
        }
    }

    std::cout << "webp_read: 成功恢复 " << diff_seq_blocks.size() << " 个序列块。" << std::endl;

    return diff_seq_blocks;
}

std::string quality_diff_base_read(const std::string& webp_base_dir, int read_length, int id) {

    // 1. 构建文件路径
    std::string parent_dir = get_webp_parent_dir(webp_base_dir.c_str());
    std::string file_path = parent_dir + "/quality_diff_base/diff_" + std::to_string(id) + ".bin";

    // 2. 打开文件
    std::ifstream ifs(file_path, std::ios::binary);
    if (!ifs.is_open()) {
        std::cerr << "diff_seq_read: 无法打开文件 " << file_path << " 进行读取。" << std::endl;
        return "";
    }

    // --- 读取文件头部 ---
    char magic[4];
    if (!ifs.read(magic, sizeof(magic)) || !(magic[0] == 'H' && magic[1] == 'u' && magic[2] == 'f' && magic[3] == 'f')) {
        std::cerr << "quality_diff_seq_read: 文件 " << file_path << " 不是有效的哈夫曼压缩文件。" << std::endl;
        ifs.close();
        return "";
    }

    uint8_t code_table_size;
    if (!ifs.read(reinterpret_cast<char*>(&code_table_size), sizeof(code_table_size))) {
        std::cerr << "quality_diff_seq_read: 读取编码表大小失败。" << std::endl;
        ifs.close();
        return "";
    }

    std::unordered_map<char, std::string> huffmanCodeTable;
    for (int i = 0; i < static_cast<int>(code_table_size); ++i) {
        char c; uint8_t code_len;
        if (!ifs.read(&c, sizeof(c)) || !ifs.read(reinterpret_cast<char*>(&code_len), sizeof(code_len))) {
            std::cerr << "quality_diff_seq_read: 读取编码表条目失败。" << std::endl;
            ifs.close();
            return "";
        }
        std::string code(code_len, '\0');
        if (!ifs.read(&code[0], code_len)) {
            std::cerr << "quality_diff_seq_read: 读取编码字符串失败。" << std::endl;
            ifs.close();
            return "";
        }
        huffmanCodeTable[c] = code;
    }

    uint64_t original_bits_length;
    if (!ifs.read(reinterpret_cast<char*>(&original_bits_length), sizeof(original_bits_length))) {
        std::cerr << "quality_diff_seq_read: 读取原始比特长度失败。" << std::endl;
        ifs.close();
        return "";
    }
    // --- 头部读取完毕 ---

    // 3. 读取压缩数据字节流
    std::vector<uint8_t> byte_data((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    if (byte_data.empty() && original_bits_length > 0) {
        std::cerr << "quality_diff_seq_read: 文件 " << file_path << " 中没有找到压缩数据。" << std::endl;
        return "";
    }

    // 4. 将字节流转换为 01 字符串
    std::string binary_string;
    binary_string.reserve(original_bits_length);
    for (uint8_t byte : byte_data) {
        for (int i = 7; i >= 0; --i) {
            binary_string += ((byte >> i) & 1) ? '1' : '0';
        }
    }

    // 5. 根据原始长度截断 01 字符串
    if (binary_string.size() > original_bits_length) {
        binary_string.resize(original_bits_length);
    }

    // 6. 重建哈夫曼树并解码
    BWT_aln::HuffmanCoder coder;
    std::shared_ptr<BWT_aln::HuffmanNode> root = coder.rebuildTreeFromCodeTable(huffmanCodeTable);
    if (!root) {
        std::cerr << "quality_diff_seq_read: 无法从编码表重建哈夫曼树。" << std::endl;
        return "";
    }
    
    // !!! 核心修改：直接返回解码后的长字符串
    std::string decoded_string = coder.decode(binary_string, root);

    if (decoded_string.empty()) {
        std::cerr << "quality_diff_seq_read: 解码失败，得到空字符串。" << std::endl;
    } else {
        std::cout << "quality_diff_seq_read: 成功解码，得到 " << decoded_string.size() << " 个字符。" << std::endl;
    }

    return decoded_string;
}

void ProcessOneID(int id, const std::string& webp_base_dir, 
                  const std::string& output_base_dir,
                  std::mutex& output_mutex,
                  std::atomic<bool>& error_occurred,
                  int read_length) {
    try {
        // 防止在已出错时继续处理
        if (error_occurred) return;
        std::vector<std::string> diff_seq_blocks = get_webp_to_diff_seq_blocks(webp_base_dir, read_length, id);
        std::string quality_diff_base = quality_diff_base_read(webp_base_dir, read_length, id);
        std::vector<std::string> quality_blocks;
        int diff_base_idx = 0;
        for(std::string& diff_seq : diff_seq_blocks)
        {
            std::string quality_score;
            quality_score.clear();
            for(auto& c:diff_seq)
            {
                if(c=='0')quality_score += 'F';
                else quality_score += quality_diff_base[diff_base_idx++];
            }
            quality_blocks.push_back(quality_score);
        }
        std::string output_file = output_base_dir + "/quality_score." + std::to_string(id) + ".bin";
        {
            std::ofstream out(output_file, std::ios::binary);
            size_t block_count = quality_blocks.size();
            out.write(reinterpret_cast<const char*>(&block_count), sizeof(size_t));
            for (const auto& block : quality_blocks) {
                size_t len = block.size();
                out.write(reinterpret_cast<const char*>(&len), sizeof(size_t));
                out.write(block.data(), len);
            }
        }
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "成功处理 ID: " << id << std::endl;
        }
    }
    catch (const std::exception& e) {
        // 设置错误标志
        error_occurred = true;
        
        // 错误输出（加锁保证线程安全）
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "处理ID " << id << " 时出错: " << e.what() << std::endl;
    }
}

// 工作线程函数
void WorkerThread(ThreadSafeQueue<int>& task_queue,
                 const std::string& webp_base_dir,
                 const std::string& output_base_dir,
                 std::mutex& output_mutex,
                 std::atomic<bool>& error_occurred,
                 int read_length) {
    int id;
    while (task_queue.Pop(id)) {
        // 如果已发生错误，提前退出
        if (error_occurred) {
            task_queue.Stop();  // 停止其他线程
            break;
        }
        ProcessOneID(id, webp_base_dir, output_base_dir, output_mutex, error_occurred, read_length);
    }
}

int webp_reconstructor_4_main(char *webp_dir, char* output_dir, int thread_n, int read_length, int total) {
        const int START_ID = 0;
        const int END_ID = 2059;  // 对应2060个文件(0-2059)
        const int NUM_THREADS = thread_n;
        
        const std::string webp_base_dir = std::string(webp_dir);
        const std::string output_base_dir = std::string(output_dir);
        
        CreateDirectoryRecursive(output_base_dir.c_str());
        
        std::cout << "使用 " << NUM_THREADS << " 个线程进行解压处理" << std::endl;
        std::cout << "处理ID范围: [" << START_ID << ", " << END_ID << "]" << std::endl;

        // 创建任务队列
        ThreadSafeQueue<int> task_queue;
        for (int id = START_ID; id <= END_ID; ++id) {
            task_queue.Push(id);
        }

        // 多线程处理控制
        std::mutex output_mutex;
        std::atomic<bool> error_occurred(false);
        std::vector<std::thread> threads;
        
        // 启动工作线程
        for (int i = 0; i < NUM_THREADS; ++i) {
            threads.emplace_back([&] {
                WorkerThread(task_queue, webp_base_dir, output_base_dir, output_mutex, error_occurred, read_length);
            });
        }

        // 等待所有线程完成
        for (auto& thread : threads) {
            thread.join();
        }

        // 检查是否发生错误
        if (error_occurred) {
            throw std::runtime_error("处理过程中出现错误，请检查日志");
        }

        std::cout << "全部文件解压验证完成" << std::endl;
        return 0;
}