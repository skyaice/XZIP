#ifndef WEBP_RECONSTRUCTOR_40_H
#define WEBP_RECONSTRUCTOR_40_H

#include <vector>
#include <string>
#include <unordered_map>
#include <atomic>
#include <mutex>
#include <chrono>
#include <queue>
#include <condition_variable>

// 全局字符与6位二进制映射表（Phred 0-39对应）
extern const std::unordered_map<char, std::string> kCharToBinary_40;
extern const std::unordered_map<std::string, char> kBinaryToChar_40;

// 进度跟踪器结构体：记录任务处理状态与时间
struct ProgressTracker_40 {
    std::atomic<int> processed{0};   // 已处理文件数
    std::atomic<int> total{0};       // 总文件数
    std::atomic<int> success{0};     // 成功处理的ID数
    std::atomic<int> errors{0};      // 失败的ID数
    std::chrono::time_point<std::chrono::steady_clock> start_time;
};

// 线程安全队列模板类（添加40后缀避免与其他头文件冲突）
template <typename T>
class ThreadSafeQueue_40 {
public:
    // 入队：添加任务到队列，唤醒等待的线程
    void Push_40(T value);

    // 出队：阻塞等待队列非空或停止标志，返回是否成功获取任务
    bool Pop_40(T& value);

    // 停止队列：设置停止标志，唤醒所有阻塞线程
    void Stop_40();

    // 检查队列是否为空（线程安全）
    bool Empty_40() const;

private:
    mutable std::mutex mutex_;               // 互斥锁：保护队列操作
    std::condition_variable cond_;           // 条件变量：实现线程阻塞/唤醒
    std::queue<T> queue_;                    // 任务队列
    bool stop_ = false;                      // 停止标志：控制队列退出
};

// -------------------------- 工具函数声明 --------------------------
/**
 * @brief 创建多级目录（支持嵌套路径）
 * @param path 目标目录路径
 */
void CreateDirectoryRecursive_40(const char* path);

/**
 * @brief 获取指定目录下所有.webp文件路径，存入vector
 * @param dir_path 目标目录路径
 * @param files 输出参数：存储.webp文件路径的vector
 * @throw std::runtime_error 目录无法打开时抛出异常
 */
void GetWebpFiles_40(const std::string& dir_path, std::vector<std::string>& files);

/**
 * @brief 自定义.webp文件排序函数（按文件名末尾数字升序）
 * @param a 第一个文件路径
 * @param b 第二个文件路径
 * @return a的数字小于b时返回true，否则返回false
 */
bool WebpFileComparator_40(const std::string& a, const std::string& b);

/**
 * @brief 从.webp文件目录恢复二维uint8矩阵（适配144序列/行、113长度/序列）
 * @param input_dir .webp文件所在目录
 * @return 恢复的二维矩阵（每行113个uint8值），目录异常时返回空矩阵
 */
std::vector<std::vector<uint8_t>> RestoreMatrixFromWebp_40(const std::string& input_dir);

/**
 * @brief 将113长度的uint8字节向量转换为Phred质量字符串
 * @param bytes 输入字节向量（长度必须为113）
 * @return 转换后的质量字符串（由kBinaryToChar_40映射生成）
 * @throw std::invalid_argument 字节向量长度非113时抛出异常
 * @throw std::runtime_error 无效二进制对时抛出异常
 */
std::string BytesToQualityString_40(const std::vector<uint8_t>& bytes);

/**
 * @brief 比较原始矩阵与恢复矩阵的一致性
 * @param original 原始矩阵
 * @param restored 恢复矩阵
 * @return 前CHECK_LIMIT行完全一致时返回true，否则返回false
 */
bool CompareMatrices_40(const std::vector<std::vector<uint8_t>>& original, 
                        const std::vector<std::vector<uint8_t>>& restored);  // 修复语法错误

/**
 * @brief 处理单个ID的完整流程：恢复矩阵、转换质量字符串、保存文件、验证矩阵
 * @param id 待处理的ID编号
 * @param webp_base_dir .webp文件基础目录（格式：xxx/）
 * @param output_base_dir 输出文件基础目录（格式：xxx/）
 * @param matrix_base_dir 原始矩阵文件基础目录
 * @param output_mutex 输出互斥锁（保证日志打印线程安全）
 * @param error_occurred 全局错误标志（原子变量，一处错误则全流程停止）
 * @param progress 进度跟踪器（更新processed/success/errors计数）
 */
void ProcessOneID_40(int id, const std::string& webp_base_dir, 
                     const std::string& output_base_dir,
                     const std::string& matrix_base_dir,
                     std::mutex& output_mutex,
                     std::atomic<bool>& error_occurred,
                     ProgressTracker_40& progress);  // 修改为带40后缀的结构体

/**
 * @brief 工作线程函数：循环从任务队列获取ID并调用ProcessOneID_40处理
 * @param task_queue 任务队列（存储待处理ID）
 * @param webp_base_dir .webp文件基础目录
 * @param output_base_dir 输出文件基础目录
 * @param matrix_base_dir 原始矩阵文件基础目录
 * @param output_mutex 输出互斥锁
 * @param error_occurred 全局错误标志
 * @param progress 进度跟踪器
 */
void WorkerThread_40(ThreadSafeQueue_40<int>& task_queue,  // 修改为带40后缀的队列
                     const std::string& webp_base_dir,
                     const std::string& output_base_dir,
                     const std::string& matrix_base_dir,
                     std::mutex& output_mutex,
                     std::atomic<bool>& error_occurred,
                     ProgressTracker_40& progress);  // 修改为带40后缀的结构体

/**
 * @brief 更新进度条显示（含百分比、处理计数、用时、剩余时间）
 * @param progress 进度跟踪器（提供processed/total/success/errors等信息）
 */
void UpdateProgressDisplay_40(const ProgressTracker_40& progress);  // 修改为带40后缀的结构体

// -------------------------- 主函数声明 --------------------------
/**
 * @brief WebP恢复器主函数：初始化任务、启动多线程、监控进度、输出报告
 * @param webp_dir .webp文件基础目录（输入，无需末尾 '/'）
 * @param output_dir 输出文件基础目录（输入，无需末尾 '/'）
 * @return 0：处理成功；1：处理失败（含异常）
 */
int webp_reconstructor_40_main(char *webp_dir, char* output_dir, int thread_n, int read_length, int total);

// -------------------------- 模板类成员函数实现 --------------------------
template <typename T>
void ThreadSafeQueue_40<T>::Push_40(T value) {  // 修正类名
    std::lock_guard<std::mutex> lock(mutex_);
    queue_.push(std::move(value));
    cond_.notify_one();
}

template <typename T>
bool ThreadSafeQueue_40<T>::Pop_40(T& value) {  // 修正类名
    std::unique_lock<std::mutex> lock(mutex_);
    cond_.wait(lock, [this]{ return !queue_.empty() || stop_; });
    
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

template <typename T>
void ThreadSafeQueue_40<T>::Stop_40() {  // 修正类名
    std::lock_guard<std::mutex> lock(mutex_);
    stop_ = true;
    cond_.notify_all();
}

template <typename T>
bool ThreadSafeQueue_40<T>::Empty_40() const {  // 修正类名
    std::lock_guard<std::mutex> lock(mutex_);
    return queue_.empty();
}

#endif // WEBP_RECONSTRUCTOR_40_H
