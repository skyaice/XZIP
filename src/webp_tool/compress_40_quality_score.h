#ifndef COMPRESS_40_QUALITY_SCORE
#define COMPRESS_40_QUALITY_SCORE

#include <string>
#include <vector>
#include <cstdint>
#include <memory>

// 单块 WebP 无损编码（块级任务单元，供线程池调用）
// matrix      : 完整质量分矩阵（只读，通过 shared_ptr 共享）
// output_dir  : 该文件的输出目录
// read_length : 每条序列的碱基数（= 图像宽度方向的像素数）
// block_idx   : 当前块编号（决定输出文件名）
// start_seq   : 该块在 matrix 中的起始行（含）
// end_seq     : 该块在 matrix 中的结束行（不含）
// image_width : 图像宽度（= sequences_per_row * read_length）
// sequences_per_row : 每行排列的序列数
void EncodeSingleWebpBlock(
    const std::vector<std::vector<uint8_t>>& matrix,
    const std::string& output_dir,
    int read_length,
    size_t block_idx,
    size_t start_seq,
    size_t end_seq,
    int image_width,
    int sequences_per_row);

// 从二进制文件读取质量分块
// id                : 文件编号
// quality_score_dir : 质量分文件所在目录
// read_length       : 每条序列的期望长度（用于校验）
std::vector<std::string> ReadQualityScoreBlocks_40(
    int id,
    char* quality_score_dir,
    int read_length);

// 主入口：多线程按块并行编码
// quality_score_dir : 质量分输入目录
// webp_dir          : WebP 输出根目录
// thread_n          : 线程池大小
// total             : 文件编号上限（处理 [0, total]，共 total+1 个文件）
// read_length       : 每条序列的碱基数
void webp_main_40(char* quality_score_dir,
    char* webp_dir,
    int thread_n,
    int total,
    int read_length,
    int webp_lossless,
    int webp_quality,
    int webp_method,
    int webp_width,
    int webp_height);

#endif  // COMPRESS_40_QUALITY_SCORE
