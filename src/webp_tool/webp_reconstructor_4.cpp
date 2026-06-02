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
#include <stdexcept>
#include <filesystem>
#include <sstream>
#include <cmath>

namespace fs = std::filesystem;

// ==================== 核心常量（与压缩侧完全对齐） ====================
const std::unordered_map<char, std::string> kCharToBinary = {
    {'F', "00"}, {':', "01"}, {'#', "10"}, {',', "11"}
};

const std::unordered_map<std::string, char> kBinaryToChar = {
    {"00", 'F'}, {"01", ':'}, {"10", '#'}, {"11", ','}
};

// ==================== 进度跟踪器（增强版） ====================
struct ProgressTracker {
    std::atomic<int> processed{0};   // 已处理文件数
    std::atomic<int> total{0};       // 总文件数
    std::atomic<int> success{0};     // 成功处理的ID数
    std::atomic<int> skipped{0};     // 跳过的ID数（目录不存在）
    std::atomic<int> errors{0};      // 失败的ID数
    std::chrono::time_point<std::chrono::steady_clock> start_time;
};

// ==================== 哈夫曼解码相关（与压缩侧对应） ====================
// 哈夫曼树节点结构（压缩侧BWT_aln::HuffmanCoder的配套解码结构）
struct HuffmanNode {
    char ch;
    uint64_t freq;
    HuffmanNode *left, *right;

    HuffmanNode(char c = '\0', uint64_t f = 0) 
        : ch(c), freq(f), left(nullptr), right(nullptr) {}

    ~HuffmanNode() {
        delete left;
        delete right;
    }
};

class HuffmanDecoder {
    public:
        HuffmanDecoder() : root_(nullptr) {}
        ~HuffmanDecoder() { delete root_; }
    
        // 从编码表重建哈夫曼树
        void buildTree(const std::unordered_map<char, std::string>& codeTable) {
            root_ = new HuffmanNode();
            for (const auto& pair : codeTable) {
                HuffmanNode* curr = root_;
                const std::string& code = pair.second;
                for (char bit : code) {
                    if (bit == '0') {
                        if (!curr->left) curr->left = new HuffmanNode();
                        curr = curr->left;
                    } else {
                        if (!curr->right) curr->right = new HuffmanNode();
                        curr = curr->right;
                    }
                }
                curr->ch = pair.first;
            }
        }
    
        // 暴露根节点（供外部逐位解码）
        HuffmanNode* root() const { return root_; }
    
    private:
        HuffmanNode* root_;
    };

// ==================== 工具函数 ====================
// 线程安全的多级目录创建（替代原CreateDirectoryRecursive）
bool CreateDirectoryRecursive(const std::string& path) {
    try {
        return fs::create_directories(path);
    } catch (const fs::filesystem_error& e) {
        std::cerr << "[Error] 创建目录失败 " << path << ": " << e.what() << std::endl;
        return false;
    }
}

// 获取目录下的webp文件列表（兼容part_*.webp/lossless_block_*.webp）
void GetWebpFiles(const std::string& dir_path, std::vector<std::string>& files) {
    files.clear();
    if (!fs::is_directory(dir_path)) {
        return; // 目录不存在直接返回空列表
    }

    for (const auto& entry : fs::directory_iterator(dir_path)) {
        if (entry.is_regular_file()) {
            const std::string filename = entry.path().filename().string();
            if (filename.size() > 5 && filename.substr(filename.size() - 5) == ".webp") {
                files.push_back(entry.path().string());
            }
        }
    }
}

// 自定义排序函数（按块编号升序）
bool WebpFileComparator(const std::string& a, const std::string& b) {
    auto extract_num = [](const std::string& path) -> int {
        std::string filename = fs::path(path).filename().string();
        size_t underscore_pos = filename.find_last_of('_');
        size_t dot_pos = filename.find_last_of('.');
        if (underscore_pos == std::string::npos || dot_pos == std::string::npos) {
            return -1;
        }
        try {
            return std::stoi(filename.substr(underscore_pos + 1, dot_pos - underscore_pos - 1));
        } catch (...) {
            return -1;
        }
    };
    return extract_num(a) < extract_num(b);
}

// 从哈夫曼编码文件恢复差异数据
std::string RestoreDiffDataFromHuffmanFile(const std::string& diff_file_path) {
    std::string diff_long_str; // 改为返回单一长字符串
    if (!fs::exists(diff_file_path)) {
        std::cerr << "[Warning] 差异数据文件不存在: " << diff_file_path << std::endl;
        return diff_long_str;
    }

    std::ifstream ifs(diff_file_path, std::ios::binary);
    if (!ifs.is_open()) {
        throw std::runtime_error("无法打开差异数据文件: " + diff_file_path);
    }

    // 1. 验证魔数
    char magic[4];
    ifs.read(magic, 4);
    if (memcmp(magic, "Huff", 4) != 0) {
        throw std::runtime_error("无效的差异数据文件格式: " + diff_file_path);
    }

    // 2. 读取编码表大小
    uint8_t code_table_size;
    ifs.read(reinterpret_cast<char*>(&code_table_size), sizeof(code_table_size));

    // 3. 读取编码表
    std::unordered_map<char, std::string> huffmanCode;
    for (int i = 0; i < code_table_size; ++i) {
        char ch;
        uint8_t code_len;
        ifs.read(&ch, sizeof(ch));
        ifs.read(reinterpret_cast<char*>(&code_len), sizeof(code_len));
        
        std::string code(code_len, '\0');
        ifs.read(&code[0], code_len);
        huffmanCode[ch] = code;
    }

    // 4. 读取原始二进制长度
    uint64_t original_bits_length;
    ifs.read(reinterpret_cast<char*>(&original_bits_length), sizeof(original_bits_length));

    // ========== 关键修正1：字节数据起始位置计算（原公式错误） ==========
    // 魔数(4) + 编码表大小(1) + 编码表条目(N*(1+1+N)) + 原始长度(8)
    // 编码表条目每个占：字符(1) + 编码长度(1) + 编码内容(N) → 动态计算
    uint64_t header_size = 4 + 1 + 8; // 魔数+表大小+原始长度
    for (const auto& pair : huffmanCode) {
        header_size += 1 + 1 + pair.second.size(); // 字符(1) + 长度(1) + 编码内容长度
    }
    // 定位到字节数据起始位置
    ifs.seekg(header_size, std::ios::beg);

    // 读取字节数据（剩余全部内容）
    std::vector<uint8_t> byte_data((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    // 5. 转换为二进制字符串
    std::string binary_str;
    binary_str.reserve(original_bits_length);
    for (uint8_t byte : byte_data) {
        for (int i = 7; i >= 0; --i) { // 高位到低位
            binary_str += (byte & (1 << i)) ? '1' : '0';
        }
    }
    binary_str.resize(original_bits_length); // 截断填充位

    // ========== 关键修正2：哈夫曼解码为单一长字符串 ==========
    HuffmanDecoder decoder;
    decoder.buildTree(huffmanCode);
    
    // 解码时不拆分序列，直接拼接成单一长字符串
    HuffmanNode* curr = decoder.root_;
    for (char bit : binary_str) {
        if (!curr) break;

        curr = (bit == '0') ? curr->left : curr->right;
        // 到达叶子节点（非结束符则拼接，结束符忽略）
        if (!curr->left && !curr->right) {
            if (curr->ch != '\0') { // 跳过结束符
                diff_long_str += curr->ch;
            }
            curr = decoder.root_; // 重置到根节点
        }
    }

    return diff_long_str;
}

// ==================== 核心解压逻辑（适配压缩侧编码规则） ====================
/**
 * 从WebP文件恢复矩阵（核心：压缩侧将F→0，其他字符→255）
 * @param input_dir WebP文件所在目录
 * @param read_length 压缩侧的read_length参数
 * @return 恢复的矩阵（0表示F，255表示其他字符）
 */
std::vector<std::vector<uint8_t>> RestoreMatrixFromWebp(const std::string& input_dir, int read_length) {
    std::vector<std::vector<uint8_t>> restored_matrix;
    if (!fs::is_directory(input_dir)) {
        return restored_matrix;
    }

    // 1. 获取并排序WebP文件
    std::vector<std::string> webp_files;
    GetWebpFiles(input_dir, webp_files);
    if (webp_files.empty()) {
        std::cerr << "[Warning] 目录 " << input_dir << " 无WebP文件" << std::endl;
        return restored_matrix;
    }
    std::sort(webp_files.begin(), webp_files.end(), WebpFileComparator);

    // 2. 计算核心参数（与压缩侧完全对齐）
    const int pixels_per_sequence = read_length;  // 每个序列占read_length个像素
    const int sequences_per_row = 16384 / pixels_per_sequence;  // 压缩侧每行序列数

    // 3. 逐文件解码
    for (const auto& file_path : webp_files) {
        // 读取WebP文件
        std::ifstream file(file_path, std::ios::binary | std::ios::ate);
        if (!file.is_open()) {
            std::cerr << "[Error] 无法打开WebP文件: " << file_path << std::endl;
            continue;
        }
        size_t file_size = file.tellg();
        file.seekg(0);
        std::vector<uint8_t> webp_data(file_size);
        if (!file.read(reinterpret_cast<char*>(webp_data.data()), file_size)) {
            std::cerr << "[Error] 读取WebP文件失败: " << file_path << std::endl;
            file.close();
            continue;
        }
        file.close();

        // WebP解码配置
        WebPDecoderConfig config;
        if (!WebPInitDecoderConfig(&config)) {
            std::cerr << "[Error] WebP解码器初始化失败: " << file_path << std::endl;
            continue;
        }

        // 获取图片信息
        VP8StatusCode status = WebPGetFeatures(webp_data.data(), webp_data.size(), &config.input);
        if (status != VP8_STATUS_OK) {
            std::cerr << "[Error] 解析WebP头失败: " << file_path << " (状态码: " << status << ")" << std::endl;
            WebPFreeDecBuffer(&config.output);
            continue;
        }

        // 解码为RGB格式（压缩侧用RGB存储）
        config.output.colorspace = MODE_RGB;
        status = WebPDecode(webp_data.data(), webp_data.size(), &config);
        if (status != VP8_STATUS_OK) {
            std::cerr << "[Error] WebP解码失败: " << file_path << " (状态码: " << status << ")" << std::endl;
            WebPFreeDecBuffer(&config.output);
            continue;
        }

        // 提取数据（压缩侧：F→0，其他→255，RGB三通道值相同）
        const int width = config.input.width;
        const int height = config.input.height;
        const uint8_t* rgb_data = config.output.u.RGB.rgb;
        const int stride = config.output.u.RGB.stride;

        // 逐行逐序列恢复
        for (int y = 0; y < height; ++y) {
            for (int seq_idx = 0; seq_idx < sequences_per_row; ++seq_idx) {
                std::vector<uint8_t> sequence(pixels_per_sequence);
                const int x_start = seq_idx * pixels_per_sequence;

                for (int px = 0; px < pixels_per_sequence; ++px) {
                    const int pixel_pos = y * stride + (x_start + px) * 3; // RGB格式：每个像素3字节
                    if (pixel_pos + 2 >= static_cast<int>(config.output.size)) {
                        sequence[px] = 0; // 边界保护
                        continue;
                    }
                    // 压缩侧R/G/B通道值相同，取R通道即可
                    sequence[px] = (rgb_data[pixel_pos] == 0) ? 0 : 255;
                }

                if (sequence.size() == pixels_per_sequence) {
                    restored_matrix.push_back(sequence);
                }
            }
        }

        WebPFreeDecBuffer(&config.output);
    }

    return restored_matrix;
}

/**
 * 从矩阵+差异数据恢复原始质量字符串
 * @param matrix 从WebP恢复的矩阵（0=F，255=其他）
 * @param diff_blocks 从哈夫曼文件恢复的差异数据
 * @param read_length 每个序列的长度
 * @return 原始质量字符串列表
 */
std::vector<std::string> RestoreQualityStrings(
    const std::vector<std::vector<uint8_t>>& matrix,
    const std::string& diff_long_str,
    int read_length) {
    
    std::vector<std::string> quality_blocks;
    if (matrix.empty()) {
        return quality_blocks;
    }

    quality_blocks.reserve(matrix.size());
    size_t diff_pos = 0; // 差异字符串的全局指针（按顺序取字符）

    for (const auto& row : matrix) {
        if (row.size() != read_length) {
            std::cerr << "[Warning] 矩阵行长度错误，预期" << read_length << "，实际" << row.size() << std::endl;
            continue;
        }

        // 初始化：所有位置默认是F
        std::string quality_str(read_length, 'F');
        
        // 遍历当前序列的每个位置，非0（非F）则从差异字符串取字符替换
        for (int i = 0; i < read_length; ++i) {
            if (row[i] != 0) { // 非F位置（压缩侧存储为255）
                if (diff_pos < diff_long_str.size()) {
                    quality_str[i] = diff_long_str[diff_pos++]; // 顺序取字符
                } else {
                    // 差异字符串长度不足，降级为F（容错）
                    std::cerr << "[Warning] 差异字符串长度不足，位置(" 
                              << quality_blocks.size() << "," << i << ") 降级为F" << std::endl;
                    quality_str[i] = 'F';
                }
            }
        }

        quality_blocks.push_back(quality_str);
    }

    // 校验：差异字符串是否完全消耗（可选，用于调试）
    if (diff_pos < diff_long_str.size()) {
        std::cerr << "[Warning] 差异字符串未完全消耗，剩余 " 
                  << diff_long_str.size() - diff_pos << " 个字符" << std::endl;
    } else if (diff_pos > diff_long_str.size()) {
        std::cerr << "[Error] 差异字符串长度不足，超出 " 
                  << diff_pos - diff_long_str.size() << " 个字符" << std::endl;
    }

    return quality_blocks;
}

// ==================== 线程安全队列（修复版） ====================
template <typename T>
class ThreadSafeQueue {
public:
    void Push(T value) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (stop_) return;
        queue_.push(std::move(value));
        cond_.notify_one();
    }

    bool Pop(T& value) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_.wait(lock, [this] { return !queue_.empty() || stop_; });
        
        if (stop_ && queue_.empty()) {
            return false;
        }
        
        if (stop_) {
            value = std::move(queue_.front());
            queue_.pop();
            return true;
        }
        
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

    size_t Size() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.size();
    }

private:
    mutable std::mutex mutex_;
    std::condition_variable cond_;
    std::queue<T> queue_;
    bool stop_ = false;
};

// ==================== 单ID处理逻辑（核心） ====================
void ProcessOneID(int id, 
                 const std::string& webp_base_dir, 
                 const std::string& output_base_dir,
                 const std::string& diff_base_dir,
                 int read_length,
                 std::mutex& output_mutex,
                 std::atomic<bool>& error_occurred,
                 ProgressTracker& progress) {
    try {
        // 错误熔断
        if (error_occurred.load()) return;

        // 1. 构建路径
        std::string webp_dir = webp_base_dir + "/combined_output." + std::to_string(id);
        std::string output_file = output_base_dir + "/quality_score." + std::to_string(id) + ".bin";
        std::string diff_file = diff_base_dir + "/diff_" + std::to_string(id) + ".bin";

        // 2. 恢复矩阵
        auto restored_matrix = RestoreMatrixFromWebp(webp_dir, read_length);
        if (restored_matrix.empty()) {
            progress.skipped++;
            progress.processed++;
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "[Skip] ID=" << id << " 无有效WebP数据" << std::endl;
            return;
        }

        // 3. 恢复差异数据
        auto diff_blocks = RestoreDiffDataFromHuffmanFile(diff_file);

        // 4. 恢复原始质量字符串
        auto quality_blocks = RestoreQualityStrings(restored_matrix, diff_blocks, read_length);
        if (quality_blocks.empty()) {
            throw std::runtime_error("恢复的质量字符串为空");
        }

        // 5. 保存为二进制文件（与压缩侧输入格式一致）
        CreateDirectoryRecursive(output_base_dir);
        std::ofstream out(output_file, std::ios::binary);
        if (!out.is_open()) {
            throw std::runtime_error("无法打开输出文件: " + output_file);
        }

        // 写入块数量 + 每个块的长度和数据
        size_t block_count = quality_blocks.size();
        out.write(reinterpret_cast<const char*>(&block_count), sizeof(size_t));
        for (const auto& block : quality_blocks) {
            size_t len = block.size();
            out.write(reinterpret_cast<const char*>(&len), sizeof(size_t));
            out.write(block.data(), len);
        }
        out.close();

        // 6. 更新进度
        progress.success++;
        progress.processed++;
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cout << "[Success] ID=" << id << " 处理完成（块数=" << quality_blocks.size() << "）" << std::endl;

    } catch (const std::exception& e) {
        error_occurred.store(true);
        progress.errors++;
        progress.processed++;
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "[Error] ID=" << id << " 处理失败: " << e.what() << std::endl;
    }
}

// ==================== 工作线程函数 ====================
void WorkerThread(ThreadSafeQueue<int>& task_queue,
                 const std::string& webp_base_dir,
                 const std::string& output_base_dir,
                 const std::string& diff_base_dir,
                 int read_length,
                 std::mutex& output_mutex,
                 std::atomic<bool>& error_occurred,
                 ProgressTracker& progress) {
    int id;
    while (task_queue.Pop(id)) {
        if (error_occurred.load()) {
            task_queue.Stop();
            break;
        }
        ProcessOneID(id, webp_base_dir, output_base_dir, diff_base_dir,
                    read_length, output_mutex, error_occurred, progress);
    }
}

// ==================== 进度显示函数 ====================
void UpdateProgressDisplay(const ProgressTracker& progress) {
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - progress.start_time);
    int processed = progress.processed.load();
    int total = progress.total.load();
    int success = progress.success.load();
    int skipped = progress.skipped.load();
    int errors = progress.errors.load();

    if (total <= 0) return;

    // 进度条
    const int bar_width = 50;
    float percentage = 100.0f * processed / total;
    int pos = static_cast<int>(bar_width * percentage / 100.0);

    std::cout << "\r[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(1) << percentage << "% ";
    std::cout << "处理: " << processed << "/" << total << " | ";
    std::cout << "成功: " << success << " | 跳过: " << skipped << " | 错误: " << errors << " | ";
    std::cout << "用时: " << elapsed.count() << "秒";

    // 剩余时间预估
    if (processed > 0 && processed < total) {
        auto eta = std::chrono::duration_cast<std::chrono::seconds>(
            elapsed * (total - processed) / processed);
        std::cout << " | 剩余: " << eta.count() << "秒";
    }

    std::cout << std::flush;
}

// ==================== 主函数（对外接口） ====================
int webp_reconstructor_4_main(char *webp_dir, char* output_dir, int thread_n, int read_length, int total) {
    try {
        // 1. 配置参数
        const int START_ID = 0;
        const int END_ID = total;
        const int NUM_THREADS = thread_n;

        // 2. 路径处理（确保末尾带/）
        std::string webp_base_dir(webp_dir);
        if (!webp_base_dir.empty() && webp_base_dir.back() != '/') {
            webp_base_dir += '/';
        }
        std::string output_base_dir(output_dir);
        if (!output_base_dir.empty() && output_base_dir.back() != '/') {
            output_base_dir += '/';
        }
        // 差异数据目录（与压缩侧一致：webp_dir的父目录/quality_diff_base）
        std::string diff_base_dir = fs::path(webp_base_dir).parent_path().string() + "/quality_diff_base/";

        // 3. 初始化进度跟踪
        ProgressTracker progress;
        progress.total = END_ID - START_ID + 1;
        progress.start_time = std::chrono::steady_clock::now();

        // 4. 打印启动信息
        std::cout << "==========================================" << std::endl;
        std::cout << "          WebP解压工具 (v4)" << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "处理范围    : ID " << START_ID << " - " << END_ID << std::endl;
        std::cout << "线程数      : " << NUM_THREADS << std::endl;
        std::cout << "Read Length : " << read_length << std::endl;
        std::cout << "源目录      : " << webp_base_dir << std::endl;
        std::cout << "输出目录    : " << output_base_dir << std::endl;
        std::cout << "差异数据目录: " << diff_base_dir << std::endl;
        std::cout << "==========================================" << std::endl;

        // 5. 创建任务队列并填充
        ThreadSafeQueue<int> task_queue;
        for (int id = START_ID; id <= END_ID; ++id) {
            task_queue.Push(id);
        }

        // 6. 启动工作线程
        std::mutex output_mutex;
        std::atomic<bool> error_occurred(false);
        std::vector<std::thread> threads;

        for (int i = 0; i < NUM_THREADS; ++i) {
            threads.emplace_back(WorkerThread,
                                std::ref(task_queue),
                                webp_base_dir,
                                output_base_dir,
                                diff_base_dir,
                                read_length,
                                std::ref(output_mutex),
                                std::ref(error_occurred),
                                std::ref(progress));
        }

        // 7. 主线程监控进度
        while (!task_queue.Empty() && !error_occurred.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            std::lock_guard<std::mutex> lock(output_mutex);
            UpdateProgressDisplay(progress);
        }

        // 8. 停止队列并等待线程结束
        task_queue.Stop();
        for (auto& t : threads) {
            if (t.joinable()) t.join();
        }

        // 9. 最终进度显示
        std::lock_guard<std::mutex> lock(output_mutex);
        UpdateProgressDisplay(progress);
        std::cout << std::endl;

        // 10. 错误检查
        if (error_occurred.load()) {
            throw std::runtime_error("处理过程中发生致命错误");
        }

        // 11. 打印总结
        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - progress.start_time);
        std::cout << "\n==========================================" << std::endl;
        std::cout << "              处理完成" << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "总任务数    : " << progress.total << std::endl;
        std::cout << "成功        : " << progress.success << std::endl;
        std::cout << "跳过        : " << progress.skipped << std::endl;
        std::cout << "错误        : " << progress.errors << std::endl;
        std::cout << "总耗时      : " << total_time.count() << " 秒" << std::endl;
        std::cout << "平均速度    : " << std::fixed << std::setprecision(2) 
                  << (progress.processed > 0 ? (double)progress.processed / total_time.count() : 0) 
                  << " ID/秒" << std::endl;
        std::cout << "==========================================" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\n[Fatal] " << e.what() << std::endl;
        return 1;
    }
}