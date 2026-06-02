// webp_converter.h
#ifndef WEBP_CONVERTER_H  // 防止头文件重复包含
#define WEBP_CONVERTER_H
#include <string>
#include <vector>
#include <cstdint>

// 声明需要被主程序调用的函数（函数原型必须和 .cpp 中的实现一致）
std::string ConvertToBinaryString(const std::string& quality_str, int read_length); 
std::vector<uint8_t> BinaryToBytes(std::string& binary, int read_lengt);
void SaveMatrixToWebpLossless(const std::vector<std::vector<uint8_t>>& matrix, const std::string& output_dir, int read_length);
void SaveAsciiMatrixToWebpLossless(const std::vector<std::vector<uint8_t>>& matrix,
    const std::string& output_path, int read_length);
void ProcessQualityBlocks(const std::vector<std::string>& quality_blocks, const std::string& output_prefix, int file_id, int read_length); 
std::vector<std::string> ReadQualityScoreBlocks(int id, char *quality_score_dir, int read_length);
void webp_main(char *quality_score_dir, char *webp_dir, int thread_n, int total, int read_length);
#endif  // WEBP_CONVERTER_H