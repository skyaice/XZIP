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

const std::unordered_map<char, std::string> kCharToBinary = {
    {'!', "000000"}, // Phred 0
    {'"', "000001"}, // Phred 1
    {'#', "000010"}, // Phred 2
    {'$', "000011"}, // Phred 3
    {'%', "000100"}, // Phred 4
    {'&', "000101"}, // Phred 5
    {'\'',"000110"}, // Phred 6 (注意单引号转义)
    {'(', "000111"}, // Phred 7
    {')', "001000"}, // Phred 8
    {'*', "001001"}, // Phred 9
    {'+', "001010"}, // Phred 10
    {',', "001011"}, // Phred 11
    {'-', "001100"}, // Phred 12
    {'.', "001101"}, // Phred 13
    {'/', "001110"}, // Phred 14
    {'0', "001111"}, // Phred 15
    {'1', "010000"}, // Phred 16
    {'2', "010001"}, // Phred 17
    {'3', "010010"}, // Phred 18
    {'4', "010011"}, // Phred 19
    {'5', "010100"}, // Phred 20
    {'6', "010101"}, // Phred 21
    {'7', "010110"}, // Phred 22
    {'8', "010111"}, // Phred 23
    {'9', "011000"}, // Phred 24
    {':', "011001"}, // Phred 25
    {';', "011010"}, // Phred 26
    {'<', "011011"}, // Phred 27
    {'=', "011100"}, // Phred 28
    {'>', "011101"}, // Phred 29
    {'?', "011110"}, // Phred 30
    {'@', "011111"}, // Phred 31
    {'A', "100000"}, // Phred 32
    {'B', "100001"}, // Phred 33
    {'C', "100010"}, // Phred 34
    {'D', "100011"}, // Phred 35
    {'E', "100100"}, // Phred 36
    {'F', "100101"}, // Phred 37
    {'G', "100110"}, // Phred 38
    {'H', "100111"}  // Phred 39
};

const std::unordered_map<std::string, char> kBinaryToChar = {
    {"000000", '!'}, // 对应 Phred 0（反转 kCharToBinary['!']）
    {"000001", '"'}, // 对应 Phred 1（反转 kCharToBinary['"']）
    {"000010", '#'}, // 对应 Phred 2（反转 kCharToBinary['#']）
    {"000011", '$'}, // 对应 Phred 3（反转 kCharToBinary['$']）
    {"000100", '%'}, // 对应 Phred 4（反转 kCharToBinary['%']）
    {"000101", '&'}, // 对应 Phred 5（反转 kCharToBinary['&']）
    {"000110", '\''},// 对应 Phred 6（反转 kCharToBinary['\'']，注意单引号转义）
    {"000111", '('}, // 对应 Phred 7（反转 kCharToBinary['(']）
    {"001000", ')'}, // 对应 Phred 8（反转 kCharToBinary[')']）
    {"001001", '*'}, // 对应 Phred 9（反转 kCharToBinary['*']）
    {"001010", '+'}, // 对应 Phred 10（反转 kCharToBinary['+']）
    {"001011", ','}, // 对应 Phred 11（反转 kCharToBinary[',']）
    {"001100", '-'}, // 对应 Phred 12（反转 kCharToBinary['-']）
    {"001101", '.'}, // 对应 Phred 13（反转 kCharToBinary['.']）
    {"001110", '/'}, // 对应 Phred 14（反转 kCharToBinary['/']）
    {"001111", '0'}, // 对应 Phred 15（反转 kCharToBinary['0']）
    {"010000", '1'}, // 对应 Phred 16（反转 kCharToBinary['1']）
    {"010001", '2'}, // 对应 Phred 17（反转 kCharToBinary['2']）
    {"010010", '3'}, // 对应 Phred 18（反转 kCharToBinary['3']）
    {"010011", '4'}, // 对应 Phred 19（反转 kCharToBinary['4']）
    {"010100", '5'}, // 对应 Phred 20（反转 kCharToBinary['5']）
    {"010101", '6'}, // 对应 Phred 21（反转 kCharToBinary['6']）
    {"010110", '7'}, // 对应 Phred 22（反转 kCharToBinary['7']）
    {"010111", '8'}, // 对应 Phred 23（反转 kCharToBinary['8']）
    {"011000", '9'}, // 对应 Phred 24（反转 kCharToBinary['9']）
    {"011001", ':'}, // 对应 Phred 25（反转 kCharToBinary[':']）
    {"011010", ';'}, // 对应 Phred 26（反转 kCharToBinary[';']）
    {"011011", '<'}, // 对应 Phred 27（反转 kCharToBinary['<']）
    {"011100", '='}, // 对应 Phred 28（反转 kCharToBinary['=']）
    {"011101", '>'}, // 对应 Phred 29（反转 kCharToBinary['>']）
    {"011110", '?'}, // 对应 Phred 30（反转 kCharToBinary['?']）
    {"011111", '@'}, // 对应 Phred 31（反转 kCharToBinary['@']）
    {"100000", 'A'}, // 对应 Phred 32（反转 kCharToBinary['A']）
    {"100001", 'B'}, // 对应 Phred 33（反转 kCharToBinary['B']）
    {"100010", 'C'}, // 对应 Phred 34（反转 kCharToBinary['C']）
    {"100011", 'D'}, // 对应 Phred 35（反转 kCharToBinary['D']）
    {"100100", 'E'}, // 对应 Phred 36（反转 kCharToBinary['E']）
    {"100101", 'F'}, // 对应 Phred 37（反转 kCharToBinary['F']）
    {"100110", 'G'}, // 对应 Phred 38（反转 kCharToBinary['G']）
    {"100111", 'H'}  // 对应 Phred 39（反转 kCharToBinary['H']）
};

// 进度跟踪器结构体
struct ProgressTracker {
    std::atomic<int> processed{0};   // 已处理文件数
    std::atomic<int> total{0};       // 总文件数
    std::atomic<int> success{0};     // 成功处理的ID数
    std::atomic<int> errors{0};      // 失败的ID数
    std::chrono::time_point<std::chrono::steady_clock> start_time;
};

// 创建多级目录
void CreateDirectoryRecursive_40(const char* path) {
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
void GetWebpFiles_40(const std::string& dir_path, std::vector<std::string>& files) {
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

std::vector<std::vector<uint8_t>> RestoreMatrixFromWebp_40(const std::string& input_dir) {
    std::vector<std::vector<uint8_t>> restored_matrix;

    // 检查目录是否存在
    struct stat info;
    if (stat(input_dir.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        std::cerr << "警告: 目录不存在，跳过处理: " << input_dir << std::endl;
        return restored_matrix; // 返回空矩阵
    }

    // 获取并排序文件列表
    std::vector<std::string> webp_files;
    try {
        GetWebpFiles_40(input_dir, webp_files);
        std::sort(webp_files.begin(), webp_files.end(), WebpFileComparator_40);
    } catch (const std::runtime_error& e) {
        std::cerr << "警告: " << e.what() << "，跳过目录: " << input_dir << std::endl;
        return restored_matrix; // 返回空矩阵
    }

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
            for (int seq_x = 0; seq_x < 144; ++seq_x) {
                std::vector<uint8_t> sequence(113);
                const int x_start = seq_x * 113;
                for (int px = 0; px < 113; ++px) {
                    const int pixel_idx = y * stride + (x_start + px) * 4;
                    sequence[px] = rgba[pixel_idx];
                }
                restored_matrix.push_back(sequence);
            }
        }
        WebPFreeDecBuffer(&config.output);
    }

    if (!restored_matrix.empty() && restored_matrix[0].size() != 113) {
        throw std::runtime_error("还原矩阵列数错误");
    }

    return restored_matrix;
}

std::string BytesToQualityString_40(const std::vector<uint8_t>& bytes) {
    if (bytes.size() != 113) {
        throw std::invalid_argument("字节向量长度必须为113");
    }

    std::string binary_str;
    binary_str.reserve(904);
    for (uint8_t byte : bytes) {
        for (int i = 7; i >= 0; --i) {
            binary_str += (byte & (1 << i)) ? '1' : '0';
        }
    }
    binary_str.resize(900); // 移除填充位

    std::string quality_str;
    for (size_t i = 0; i < 900; i += 6) {
        std::string pair = binary_str.substr(i, 6);
        auto it = kBinaryToChar.find(pair);
        if (it == kBinaryToChar.end()) {
            throw std::runtime_error("无效的二进制对: " + pair);
        }
        quality_str += it->second;
    }

    return quality_str;
}

bool CompareMatrices_40(const std::vector<std::vector<uint8_t>>& original, const std::vector<std::vector<uint8_t>>& restored) {
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
    void Push_40(T value) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push(std::move(value));
        cond_.notify_one();
    }

    bool Pop_40(T& value) {
        std::unique_lock<std::mutex> lock(mutex_);
        // 等待条件：队列非空 或 停止标志为真
        cond_.wait(lock, [this]{ return !queue_.empty() || stop_; });
        if (stop_ && queue_.empty()) { // 新增：队列空且停止，才返回false
            return false;
        }
        if (stop_) { // 停止标志为真但队列非空，仍处理剩余任务
            value = std::move(queue_.front());
            queue_.pop();
            return true;
        }
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
                  const std::string& matrix_base_dir,
                  std::mutex& output_mutex,
                  std::atomic<bool>& error_occurred,
                  ProgressTracker& progress) {
    try {
        // 防止在已出错时继续处理
        if (error_occurred.load()) return;

        // 恢复矩阵
        std::string webp_dir = webp_base_dir + "combined_output." + std::to_string(id);
        auto restored_matrix = RestoreMatrixFromWebp_40(webp_dir);

        // 如果矩阵为空（目录不存在），跳过该ID的处理
        if (restored_matrix.empty()) {
            progress.errors++;
            progress.processed++;
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "跳过缺失目录: ID " << id << std::endl;
            return;
        }

        // 转换质量块
        std::vector<std::string> quality_blocks;
        for (const auto& row : restored_matrix) {
            quality_blocks.push_back(BytesToQualityString_40(row));
        }

        // 保存质量文件
        std::string output_file = output_base_dir + "quality_score." + std::to_string(id) + ".bin";
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

        // 验证矩阵
        std::string matrix_file = matrix_base_dir + "matrix_" + std::to_string(id) + ".bin";
        std::ifstream in(matrix_file, std::ios::binary);
        uint32_t rows, cols;
        in.read(reinterpret_cast<char*>(&rows), sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(&cols), sizeof(uint32_t));
        std::vector<std::vector<uint8_t>> original_matrix(rows, std::vector<uint8_t>(cols));
        for (auto& row : original_matrix) {
            in.read(reinterpret_cast<char*>(row.data()), cols);
        }

        if (!CompareMatrices_40(original_matrix, restored_matrix)) {
            throw std::runtime_error("矩阵验证失败");
        }

        // 更新成功计数
        progress.success++;
        
        // 输出结果（加锁保证线程安全）
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            // 不再打印每个ID的成功信息，避免干扰进度条
        }
    }
    catch (const std::exception& e) {
        // 设置错误标志
        error_occurred = true;
        
        // 更新错误计数
        progress.errors++;
        
        // 错误输出（加锁保证线程安全）
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "\n处理ID " << id << " 时出错: " << e.what() << std::endl;
    }
    
    // 更新已处理计数
    progress.processed++;
}

// 工作线程函数
void WorkerThread_40(ThreadSafeQueue<int>& task_queue,
                 const std::string& webp_base_dir,
                 const std::string& output_base_dir,
                 const std::string& matrix_base_dir,
                 std::mutex& output_mutex,
                 std::atomic<bool>& error_occurred,
                 ProgressTracker& progress) {
    int id;
    while (task_queue.Pop_40(id)) {
        // 如果已发生错误，提前退出
        if (error_occurred.load()) {
            task_queue.Stop_40();  // 停止其他线程
            break;
        }
        ProcessOneID_40(id, webp_base_dir, output_base_dir, matrix_base_dir, output_mutex, error_occurred, progress);
    }
}

// 更新进度显示
void UpdateProgressDisplay_40(const ProgressTracker& progress) {
    auto now = std::chrono::steady_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(now - progress.start_time);
    auto elapsed_millis = std::chrono::duration_cast<std::chrono::milliseconds>(now - progress.start_time);
    
    int processed = progress.processed.load();
    int total = progress.total.load();
    int success = progress.success.load();
    int errors = progress.errors.load();
    
    if (total <= 0) return;
    
    // 计算进度百分比
    float percentage = 100.0f * processed / total;
    
    // 进度条显示
    const int bar_width = 50;
    int pos = static_cast<int>(bar_width * percentage / 100.0);
    
    std::cout << "\r[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] ";
    
    // 数字统计
    std::cout << std::setw(5) << std::fixed << std::setprecision(1) << percentage << "% "
              << "处理: " << processed << "/" << total << " "
              << "成功: " << success << " "
              << "错误: " << errors << " "
              << "用时: " << elapsed_seconds.count() << "秒";
    
    // 计算剩余时间（如果可能）
    if (processed > 0 && processed < total) {
        auto remaining = std::chrono::duration_cast<std::chrono::seconds>(
            (now - progress.start_time) * (total - processed) / processed);
        std::cout << " 剩余: " << remaining.count() << "秒";
    }
    
    // 刷新输出
    std::cout << std::flush;
}

int webp_reconstructor_40_main(char *webp_dir, char* output_dir) {
    try {
        const int START_ID = 0;
        const int END_ID = 2059; // 对应2060个文件(0-2059)
        const int NUM_THREADS = 36;
        
        // 确保目录路径正确格式
        std::string webp_base_dir = webp_dir;
        if (!webp_base_dir.empty() && webp_base_dir.back() != '/') {
            webp_base_dir += '/';
        }
        
        std::string output_base_dir = output_dir;
        if (!output_base_dir.empty() && output_base_dir.back() != '/') {
            output_base_dir += '/';
        }
        
        // 验证用路径（根据实际情况调整）
        const std::string matrix_base_dir = "/home/user/yexiang/save_image_matrix/";
        
        // 创建输出目录
        CreateDirectoryRecursive_40(output_base_dir.c_str());
        
        // 初始化进度跟踪器
        ProgressTracker progress;
        progress.total = END_ID - START_ID + 1;
        progress.start_time = std::chrono::steady_clock::now();
        
        // 打印启动信息
        std::time_t start_time_t = std::chrono::system_clock::to_time_t(
            std::chrono::system_clock::now());
        std::cout << "使用 " << NUM_THREADS << " 个线程进行解压处理" << std::endl;
        std::cout << "处理ID范围: [" << START_ID << ", " << END_ID << "]" << std::endl;
        std::cout << "开始时间: " << std::ctime(&start_time_t);
        std::cout << "源目录: " << webp_base_dir << std::endl;
        std::cout << "输出目录: " << output_base_dir << std::endl;

        // 创建任务队列
        ThreadSafeQueue<int> task_queue;
        for (int id = START_ID; id <= END_ID; ++id) {
            task_queue.Push_40(id);
        }

        test_queue.Stop();

        // 多线程处理控制
        std::mutex output_mutex;
        std::atomic<bool> error_occurred(false);
        std::vector<std::thread> threads;
        
        // 启动工作线程
        for (int i = 0; i < NUM_THREADS; ++i) {
            threads.emplace_back([&] {
                WorkerThread_40(task_queue, webp_base_dir, output_base_dir, matrix_base_dir, 
                            output_mutex, error_occurred, progress);
            });
        }

        // 主线程监控进度
        while (!task_queue.Empty_40() && !error_occurred.load_40()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            std::lock_guard<std::mutex> lock(output_mutex);
            UpdateProgressDisplay_40(progress);
        }
        
        // 等待所有线程完成
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
        
        // 更新最终进度
        UpdateProgressDisplay_40(progress);
        std::cout << std::endl;

        // 检查是否发生错误
        if (error_occurred.load()) {
            throw std::runtime_error("处理过程中出现错误");
        }

        // 打印最终报告
        auto end_time = std::chrono::steady_clock::now();
        auto total_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - progress.start_time);
        auto total_millis = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - progress.start_time);
        std::time_t end_time_t = std::chrono::system_clock::to_time_t(
            std::chrono::system_clock::now());
        
        std::cout << "\n==========================================" << std::endl;
        std::cout << "        所有文件处理完成" << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << "总任务数: " << progress.total << std::endl;
        std::cout << "成功处理: " << progress.success << std::endl;
        std::cout << "失败任务: " << progress.errors << std::endl;
        std::cout << "开始时间: " << std::ctime(&start_time_t);
        std::cout << "结束时间: " << std::ctime(&end_time_t);
        std::cout << "总耗时  : " << total_seconds.count() << " 秒 (" 
                  << total_millis.count() << " 毫秒)" << std::endl;
        std::cout << "处理速度: " << std::fixed << std::setprecision(2) 
                  << (progress.total * 1000.0 / total_millis.count()) << " 任务/秒" << std::endl;
        std::cout << "==========================================" << std::endl;
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\n程序错误: " << e.what() << std::endl;
        return 1;
    }
}