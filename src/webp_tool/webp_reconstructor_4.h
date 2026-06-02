#ifndef WEBP_RECONSTRUCTOR_4_H
#define WEBP_RECONSTRUCTOR_4_H

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


// --------------------------- 全局常量声明 ---------------------------
// 字符与二进制的映射表（外部引用声明，定义在 .cpp 中）
extern const std::unordered_map<char, std::string> kCharToBinary;
extern const std::unordered_map<std::string, char> kBinaryToChar;


// --------------------------- 线程安全队列类声明 ---------------------------
// 模板类：用于多线程任务调度的线程安全队列
template <typename T>
class ThreadSafeQueue {
public:
    // 向队列推送任务（移动语义，避免拷贝）
    void Push(T value);

    // 从队列弹出任务（阻塞直到有任务或停止信号）
    // 返回值：true 表示成功弹出任务，false 表示队列已停止且无任务
    bool Pop(T& value);

    // 停止队列（唤醒所有阻塞的线程，后续 Pop 会返回 false）
    void Stop();

    // 检查队列是否为空（线程安全）
    bool Empty() const;

private:
    mutable std::mutex mutex_;          // 互斥锁（mutable 允许 const 成员函数修改）
    std::condition_variable cond_;      // 条件变量（用于线程唤醒）
    std::queue<T> queue_;               // 底层任务队列
    bool stop_ = false;                 // 停止标志（控制队列生命周期）
};


// --------------------------- 工具函数声明 ---------------------------
/**
 * @brief 创建多级目录（递归创建，支持 Linux/macOS 路径格式）
 * @param path 目标目录路径（C 风格字符串）
 */
void CreateDirectoryRecursive(const char* path);

/**
 * @brief 获取指定目录下所有 .webp 文件的路径列表
 * @param dir_path 目标目录路径
 * @param files 输出参数：存储 .webp 文件路径的向量
 * @throw std::runtime_error 目录无法打开时抛出异常
 */
void GetWebpFiles(const std::string& dir_path, std::vector<std::string>& files);

/**
 * @brief WebP 文件排序比较器（按文件名中最后一个 "_" 后的数字升序排序）
 * @param a 第一个文件路径
 * @param b 第二个文件路径
 * @return true 表示 a 应排在 b 前面，false 反之
 */
bool WebpFileComparator(const std::string& a, const std::string& b);

/**
 * @brief 从指定目录的 WebP 文件中恢复图像矩阵
 * @param input_dir 包含 WebP 文件的目录路径
 * @return 恢复后的图像矩阵（每行是一个 uint8_t 向量，每行长度固定为 38）
 * @throw std::runtime_error 解码失败、文件读取错误或矩阵格式错误时抛出异常
 */
std::vector<std::vector<uint8_t>> RestoreMatrixFromWebpTiled(const std::string& file_path, int pixels_per_sequence);
std::vector<std::vector<uint8_t>> RestoreMatrixFromWebp(const std::string& input_dir);

/**
 * @brief 将 38 字节的向量转换为质量字符串（150 个字符，由 F/:/, #/, , 组成）
 * @param bytes 输入字节向量（长度必须为 38）
 * @return 转换后的质量字符串
 * @throw std::invalid_argument 输入字节向量长度不为 38 时抛出异常
 * @throw std::runtime_error 二进制对无法匹配时抛出异常
 */
std::string BytesToQualityString(const std::vector<uint8_t>& bytes);

/**
 * @brief 比较两个矩阵的一致性（前 N 行，N 为原始矩阵行数）
 * @param original 原始矩阵（用于验证的基准）
 * @param restored 恢复后的矩阵（待验证的矩阵）
 * @return true 表示前 N 行完全一致，false 反之（差异信息会输出到 stderr）
 */
bool CompareMatrices(const std::vector<std::vector<uint8_t>>& original, const std::vector<std::vector<uint8_t>>& restored);


// --------------------------- 任务处理函数声明 ---------------------------
/**
 * @brief 处理单个 ID 的任务（恢复矩阵、转换质量块、保存文件、验证一致性）
 * @param id 任务 ID（对应 combined_output.{id} 目录和 quality_score.{id}.bin 文件）
 * @param webp_base_dir WebP 文件基础目录（combined_output.{id} 的父目录）
 * @param output_base_dir 输出文件基础目录（quality_score.{id}.bin 的父目录）
 * @param matrix_base_dir 原始矩阵文件基础目录（matrix_{id}.bin 的父目录）
 * @param output_mutex 输出互斥锁（保证多线程日志输出线程安全）
 * @param error_occurred 错误标志（原子变量，任一任务出错时设为 true）
 */
void ProcessOneID(int id, const std::string& webp_base_dir, 
                  const std::string& output_base_dir,
                  const std::string& matrix_base_dir,
                  std::mutex& output_mutex,
                  std::atomic<bool>& error_occurred);

/**
 * @brief 工作线程函数（从任务队列获取 ID 并调用 ProcessOneID 处理）
 * @param task_queue 线程安全任务队列（存储待处理的 ID）
 * @param webp_base_dir WebP 文件基础目录
 * @param output_base_dir 输出文件基础目录
 * @param matrix_base_dir 原始矩阵文件基础目录
 * @param output_mutex 输出互斥锁
 * @param error_occurred 错误标志（任一任务出错时停止所有线程）
 */
void WorkerThread(ThreadSafeQueue<int>& task_queue,
                 const std::string& webp_base_dir,
                 const std::string& output_base_dir,
                 const std::string& matrix_base_dir,
                 std::mutex& output_mutex,
                 std::atomic<bool>& error_occurred);



int webp_reconstructor_4_main(char *webp_dir, char* output_dir, int thread_n, int read_length, int total);


#endif  // WEBP_RECONSTRUCTOR_4_H