/*
 * BWT_aln.cpp
 *
 *  Created on: 2022年5月26日
 *      Author: fenghe
 */

 #include <string.h>
 #include <getopt.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <libgen.h>
 #include <sys/time.h>
 #include <sys/stat.h>
 #include <sys/types.h>
 #include <iostream>
 #include <sstream>
 #include <queue>
 #include <algorithm>
 #include <inttypes.h>
 #include <filesystem>
 #include <system_error>
 #include <chrono>
 #include <ctime>
 #include <random>
 #include <array>
 #include <fstream>
 #include <zstd.h>
 #include "BWT_aln.hpp"
 #include "CPPLIB/tools.hpp"
 extern "C"{
 #include "clib/kthread.h"
 #include "clib/utils.h"
 #include "clib/desc.h"
 }
 #include "occ_RST_DEF.hpp"
 #include "./BWT_idx/var_map.hpp"
 #include "./webp_tool/webp_converter.h"
 #include "./webp_tool/compress_40_quality_score.h"
 #include "./webp_tool/webp_reconstructor_4.h"
 #include "./webp_tool/webp_reconstructor_40.h"
 
 
 #include <set>
 
 #include <math.h>
 #include <zlib.h>
 #include "./htslib/htslib/sam.h"

 #ifndef HTS_LINE_INCLUDE_NL
 #define HTS_LINE_INCLUDE_NL 1  // 旧版本中用1表示包含换行符
 #endif

 #ifndef KS_SEP_LINE
 #define KS_SEP_LINE 2  // 旧版本中KS_SEP_LINE的典型值
 #endif 
 
 using namespace BWT_aln;
 using FastaData = std::vector<std::pair<std::string, std::string>>;
 namespace fs = std::filesystem;

 using PerfClock = std::chrono::steady_clock;

inline double elapsed_ms(const PerfClock::time_point& start,
                         const PerfClock::time_point& end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

struct ProducerPerf {
    size_t batch_cnt = 0;
    size_t read_pairs = 0;
    double batch_fill_ms = 0.0;   // 从 batch 开始装到装满
    double push_ms = 0.0;         // push_batch 阻塞/入队耗时
};

struct ConsumerPerf {
    size_t batch_cnt = 0;
    size_t pair_cnt = 0;
    size_t error_pair_cnt = 0;
    size_t flush_cnt = 0;

    double pop_wait_ms = 0.0;     // 等队列拿到 batch 的时间
    double decode_ms = 0.0;       // get_bam_seq + qual_str
    double ref_ms = 0.0;          // window_id/ref 构造
    double error_io_ms = 0.0;     // write_error_reads_to_fq
    double hamming_ms = 0.0;      // min_hamming_distance
    double hap_ms = 0.0;          // find_best_hap
    double block_ms = 0.0;        // 填充 Compress_block + push local_blocks
    double flush_ms = 0.0;        // flush_local_blocks
    double cleanup_ms = 0.0;      // bam_destroy1 等
    double total_batch_ms = 0.0;  // 整个 batch 总耗时
};

static std::mutex g_perf_print_mutex;
 

 #define PIPELINE_T_NUM 3//one for reading; one for classifying; one for writing
 #define STEP_NUM PIPELINE_T_NUM
 #define N_NEEDED_B 2000000 //2M read pair per time
 #define N_NEEDED N_NEEDED_B+100
 #define WINDOW_ID_MAX 50000000
 #define MAX_MEM_BUFFER_SIZE 1000000
uint64_t read_sum = 0;
std::atomic<uint64_t> g_error_qname_counter{0};
int cnt_14;
 int g_global_read_length = 0;
 std::vector<int>chr_bg_wb_ID;
 std::vector<Window_t> g_window_info(WINDOW_ID_MAX);
 //debug below
 int max_match_length_count[455];
 int hamming_count[455];
 int hamming_origin[455];
 int diff_base_length[4000][455];
 uint64_t quality_score_count[256];
 int invalid_window = 0;
 //debug above
 
 FastaData genome;

 char charTchar[]={
	 /*   0 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /*  16 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /*  32 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /*  48 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /*  64 */ 'T', 'A', 'G', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /*    	   A         C                   G                                  N */
	 /*  80 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',//'Z'
	 /*                        T */
	 /*  96 */ 'T', 'A', 'G', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /*         a         c                   g                                  n */
	 /* 112 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /*                        t */
	 /* 128 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /* 144 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /* 160 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /* 176 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /* 192 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /* 208 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /* 224 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
	 /* 240 */ 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
 };
 
	 static unsigned char seq_nt4_table[256] = {
		 0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		 4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	 };
 
// 在原有CharPtrHash结构体后添加std::string重载
struct CharPtrHash {
    // 原有：支持const char*
    size_t operator()(const char* str) const {
        size_t hash = 0;
        while (*str) {
            hash = hash * 31 + *str++;
        }
        return hash;
    }
    // 新增：支持std::string（调用原有const char*版本）
    size_t operator()(const std::string& str) const {
        return operator()(str.c_str()); // 复用原有哈希逻辑
    }
};

// 在原有CharPtrEqual结构体后添加std::string重载
struct CharPtrEqual {
    // 原有：支持const char*
    bool operator()(const char* a, const char* b) const {
        return strcmp(a, b) == 0;
    }
    // 新增：支持std::string（调用原有const char*版本）
    bool operator()(const std::string& a, const std::string& b) const {
        return operator()(a.c_str(), b.c_str()); // 复用原有比较逻辑
    }
};
 
 
 size_t bitLength(uint64_t num) {
	 size_t length = 0;
	 while (num > 0) {
		 length += 2;
		 num >>= 2; // 每次右移两位
	 }
	 return length;
 }
 
 std::string convertBitsToString(const std::vector<unsigned char>& binary_data, size_t code_len) {
    std::string bit_str;
    bit_str.reserve(code_len);  // 预分配空间
    
    for (size_t i = 0; i < code_len; ++i) {
        size_t byte_idx = i / 8;
        int bit_offset = 7 - (i % 8);  // 保持与写入时一致的高位优先
        
        // 检查是否超出binary_data范围（防止越界）
        if (byte_idx >= binary_data.size()) {
            break;
        }
        
        // 提取对应位的值
        unsigned char byte = binary_data[byte_idx];
        bool is_set = (byte & (1 << bit_offset)) != 0;
        bit_str += is_set ? '1' : '0';
    }
    
    return bit_str;
}

bool execute_system_command(const std::string& cmd, const std::string& err_msg) {
    fprintf(stderr, "Executing command: %s\n", cmd.c_str());
    int ret = system(cmd.c_str());
    if (ret != 0) {
        fprintf(stderr, "ERROR: %s (return code: %d)\n", err_msg.c_str(), ret);
        return false;
    }
    return true;
}

static std::string shell_quote(const std::string& value) {
    std::string quoted = "'";
    for (char c : value) {
        if (c == '\'') quoted += "'\\''";
        else quoted += c;
    }
    quoted += "'";
    return quoted;
}

bool directory_exists(const std::string& dir_path) {
    std::error_code ec;
    return fs::is_directory(dir_path, ec);
}

bool file_exists(const std::string& file_path) {
    std::error_code ec;
    return fs::is_regular_file(file_path, ec);
}

template <typename T>
static bool write_binary_value(std::ostream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return out.good();
}

template <typename T>
static bool read_binary_value(std::istream& in, T& value) {
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
    return in.good();
}

static bool ensure_parent_directory(const fs::path& file_path) {
    std::error_code ec;
    fs::path parent = file_path.parent_path();
    if (parent.empty()) return true;
    fs::create_directories(parent, ec);
    return !ec;
}

static bool zstd_compress_stream_to_stream(const fs::path& input_path,
                                           std::ostream& out,
                                           int thread_n,
                                           int zstd_level,
                                           uint64_t& compressed_size) {
    compressed_size = 0;
    std::ifstream in(input_path, std::ios::binary);
    if (!in) {
        std::cerr << "[XZIP] cannot open input file for zstd: "
                  << input_path << std::endl;
        return false;
    }

    ZSTD_CCtx* cctx = ZSTD_createCCtx();
    if (cctx == nullptr) {
        std::cerr << "[XZIP] ZSTD_createCCtx failed" << std::endl;
        return false;
    }

    const int level = std::clamp(zstd_level, 1, 19);
    size_t zr = ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, level);
    if (ZSTD_isError(zr)) {
        std::cerr << "[XZIP] zstd level error: " << ZSTD_getErrorName(zr) << std::endl;
        ZSTD_freeCCtx(cctx);
        return false;
    }
    zr = ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, std::max(0, thread_n));
    if (ZSTD_isError(zr)) {
        std::cerr << "[XZIP] zstd worker error: " << ZSTD_getErrorName(zr) << std::endl;
        ZSTD_freeCCtx(cctx);
        return false;
    }

    std::vector<char> in_buf(ZSTD_CStreamInSize());
    std::vector<char> out_buf(ZSTD_CStreamOutSize());

    bool last_chunk = false;
    while (!last_chunk) {
        in.read(in_buf.data(), static_cast<std::streamsize>(in_buf.size()));
        size_t read_size = static_cast<size_t>(in.gcount());
        last_chunk = in.eof();

        ZSTD_inBuffer input = { in_buf.data(), read_size, 0 };
        ZSTD_EndDirective mode = last_chunk ? ZSTD_e_end : ZSTD_e_continue;
        bool finished = false;
        while (!finished) {
            ZSTD_outBuffer output = { out_buf.data(), out_buf.size(), 0 };
            size_t remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
            if (ZSTD_isError(remaining)) {
                std::cerr << "[XZIP] zstd compress error: "
                          << ZSTD_getErrorName(remaining) << std::endl;
                ZSTD_freeCCtx(cctx);
                return false;
            }
            if (output.pos > 0) {
                out.write(out_buf.data(), static_cast<std::streamsize>(output.pos));
                compressed_size += static_cast<uint64_t>(output.pos);
                if (!out.good()) {
                    ZSTD_freeCCtx(cctx);
                    return false;
                }
            }
            finished = last_chunk ? (remaining == 0) : (input.pos == input.size);
        }
    }

    ZSTD_freeCCtx(cctx);
    return true;
}

static bool zstd_decompress_stream_from_stream(std::istream& in,
                                               uint64_t compressed_size,
                                               const fs::path& output_path) {
    if (!ensure_parent_directory(output_path)) {
        std::cerr << "[XZIP] cannot create output parent for "
                  << output_path << std::endl;
        return false;
    }
    std::ofstream out(output_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        std::cerr << "[XZIP] cannot open output file: " << output_path << std::endl;
        return false;
    }

    ZSTD_DCtx* dctx = ZSTD_createDCtx();
    if (dctx == nullptr) {
        std::cerr << "[XZIP] ZSTD_createDCtx failed" << std::endl;
        return false;
    }

    std::vector<char> in_buf(ZSTD_DStreamInSize());
    std::vector<char> out_buf(ZSTD_DStreamOutSize());
    uint64_t remain = compressed_size;

    while (remain > 0) {
        size_t to_read = static_cast<size_t>(
            std::min<uint64_t>(remain, static_cast<uint64_t>(in_buf.size())));
        in.read(in_buf.data(), static_cast<std::streamsize>(to_read));
        size_t read_size = static_cast<size_t>(in.gcount());
        if (read_size == 0) {
            std::cerr << "[XZIP] truncated zstd stream while reading "
                      << output_path << std::endl;
            ZSTD_freeDCtx(dctx);
            return false;
        }
        remain -= read_size;

        ZSTD_inBuffer input = { in_buf.data(), read_size, 0 };
        while (input.pos < input.size) {
            ZSTD_outBuffer output = { out_buf.data(), out_buf.size(), 0 };
            size_t zr = ZSTD_decompressStream(dctx, &output, &input);
            if (ZSTD_isError(zr)) {
                std::cerr << "[XZIP] zstd decompress error: "
                          << ZSTD_getErrorName(zr) << std::endl;
                ZSTD_freeDCtx(dctx);
                return false;
            }
            if (output.pos > 0) {
                out.write(out_buf.data(), static_cast<std::streamsize>(output.pos));
                if (!out.good()) {
                    ZSTD_freeDCtx(dctx);
                    return false;
                }
            }
        }
    }

    ZSTD_freeDCtx(dctx);
    return true;
}

static bool write_xzip_entry(std::ofstream& archive,
                             const fs::path& base_dir,
                             const fs::path& file_path,
                             int thread_n,
                             int zstd_level) {
    std::error_code ec;
    if (!fs::is_regular_file(file_path, ec)) return true;

    std::string rel_path = fs::relative(file_path, base_dir, ec).generic_string();
    if (ec || rel_path.empty()) {
        std::cerr << "[XZIP] cannot make relative path for " << file_path << std::endl;
        return false;
    }
    if (rel_path.size() > UINT16_MAX) {
        std::cerr << "[XZIP] archive path too long: " << rel_path << std::endl;
        return false;
    }

    uint16_t path_len = static_cast<uint16_t>(rel_path.size());
    uint64_t original_size = static_cast<uint64_t>(fs::file_size(file_path, ec));
    if (ec) {
        std::cerr << "[XZIP] cannot stat " << file_path << std::endl;
        return false;
    }

    if (!write_binary_value(archive, path_len)) return false;
    archive.write(rel_path.data(), path_len);
    if (!write_binary_value(archive, original_size)) return false;

    std::streampos compressed_size_pos = archive.tellp();
    uint64_t compressed_size = 0;
    if (!write_binary_value(archive, compressed_size)) return false;

    if (!zstd_compress_stream_to_stream(file_path, archive, thread_n,
                                        zstd_level, compressed_size)) {
        return false;
    }

    std::streampos end_pos = archive.tellp();
    archive.seekp(compressed_size_pos);
    if (!write_binary_value(archive, compressed_size)) return false;
    archive.seekp(end_pos);

    std::cout << "[XZIP] packed " << rel_path
              << " raw=" << original_size
              << " zstd=" << compressed_size << std::endl;
    return archive.good();
}

static bool create_xzip_archive(const std::string& parent_dir,
                                const std::string& archive_path,
                                const std::vector<std::string>& paths,
                                int thread_n,
                                int zstd_level) {
    fs::path base(parent_dir);
    std::ofstream archive(archive_path, std::ios::binary | std::ios::trunc);
    if (!archive) {
        std::cerr << "[XZIP] cannot create archive: " << archive_path << std::endl;
        return false;
    }

    const char magic[8] = { 'X', 'Z', 'I', 'P', 'A', 'R', '1', '\0' };
    archive.write(magic, sizeof(magic));
    uint32_t version = 1;
    if (!write_binary_value(archive, version)) return false;

    uint32_t entry_count = 0;
    std::streampos count_pos = archive.tellp();
    if (!write_binary_value(archive, entry_count)) return false;

    for (const std::string& rel : paths) {
        fs::path path = base / rel;
        std::error_code ec;
        if (fs::is_regular_file(path, ec)) {
            if (!write_xzip_entry(archive, base, path, thread_n, zstd_level)) return false;
            entry_count++;
        } else if (fs::is_directory(path, ec)) {
            for (const auto& item : fs::recursive_directory_iterator(path, ec)) {
                if (ec) break;
                if (item.is_regular_file(ec)) {
                    if (!write_xzip_entry(archive, base, item.path(), thread_n, zstd_level)) return false;
                    entry_count++;
                }
            }
            if (ec) {
                std::cerr << "[XZIP] cannot iterate directory: " << path << std::endl;
                return false;
            }
        } else {
            std::cerr << "[XZIP] skip missing path: " << path << std::endl;
        }
    }

    archive.seekp(count_pos);
    if (!write_binary_value(archive, entry_count)) return false;
    archive.seekp(0, std::ios::end);
    archive.close();

    std::cout << "[XZIP] archive written: " << archive_path
              << " entries=" << entry_count << std::endl;
    return true;
}

static bool restore_xzip_archive(const std::string& parent_dir,
                                 const std::string& archive_path) {
    std::ifstream archive(archive_path, std::ios::binary);
    if (!archive) {
        std::cerr << "[XZIP] archive not found: " << archive_path << std::endl;
        return false;
    }

    char magic[8] = {};
    archive.read(magic, sizeof(magic));
    if (memcmp(magic, "XZIPAR1", 7) != 0) {
        std::cerr << "[XZIP] invalid archive magic: " << archive_path << std::endl;
        return false;
    }

    uint32_t version = 0;
    uint32_t entry_count = 0;
    if (!read_binary_value(archive, version) ||
        !read_binary_value(archive, entry_count) ||
        version != 1) {
        std::cerr << "[XZIP] invalid archive header: " << archive_path << std::endl;
        return false;
    }

    fs::path base(parent_dir);
    for (uint32_t i = 0; i < entry_count; ++i) {
        uint16_t path_len = 0;
        uint64_t original_size = 0;
        uint64_t compressed_size = 0;
        if (!read_binary_value(archive, path_len)) return false;
        std::string rel_path(path_len, '\0');
        archive.read(rel_path.data(), path_len);
        if (!archive.good() ||
            !read_binary_value(archive, original_size) ||
            !read_binary_value(archive, compressed_size)) {
            return false;
        }

        fs::path output_path = base / fs::path(rel_path);
        if (!zstd_decompress_stream_from_stream(archive, compressed_size, output_path)) {
            return false;
        }

        std::error_code ec;
        uint64_t restored_size = static_cast<uint64_t>(fs::file_size(output_path, ec));
        if (ec || restored_size != original_size) {
            std::cerr << "[XZIP] restored size mismatch for " << rel_path
                      << " expected=" << original_size
                      << " got=" << restored_size << std::endl;
            return false;
        }
    }

    std::cout << "[XZIP] archive restored: " << archive_path
              << " entries=" << entry_count << std::endl;
    return true;
}

bool tar_zstd_directory(const std::string& parent_dir,
                        const std::string& dir_name,
                        int thread_n,
                        int zstd_level) {
    const std::string tar_path = parent_dir + "/" + dir_name + ".tar";
    const int zstd_threads = std::max(1, thread_n);
    const int level = std::clamp(zstd_level, 1, 19);

    std::string tar_cmd = "tar -C " + shell_quote(parent_dir) +
                          " -cf " + shell_quote(tar_path) +
                          " " + shell_quote(dir_name);
    if (!execute_system_command(tar_cmd, "tar package failed: " + dir_name)) {
        return false;
    }

    std::string zstd_cmd = "zstd -T" + std::to_string(zstd_threads) +
                           " -" + std::to_string(level) +
                           " -f " + shell_quote(tar_path);
    return execute_system_command(zstd_cmd, "zstd compress failed: " + tar_path);
}

bool zstd_compress_file(const std::string& file_path,
                        int thread_n,
                        int zstd_level,
                        const std::string& err_msg) {
    std::ofstream out(file_path + ".zst", std::ios::binary | std::ios::trunc);
    if (!out) {
        std::cerr << "ERROR: cannot create " << file_path << ".zst" << std::endl;
        return false;
    }
    uint64_t compressed_size = 0;
    bool ok = zstd_compress_stream_to_stream(file_path, out, thread_n,
                                             zstd_level, compressed_size);
    if (!ok) {
        std::cerr << "ERROR: " << err_msg << std::endl;
        return false;
    }
    std::cout << "[ZSTD] compressed " << file_path
              << " -> " << file_path << ".zst"
              << " bytes=" << compressed_size << std::endl;
    return true;
}

bool restore_zstd_tar_directory(const std::string& parent_dir,
                                const std::string& dir_name) {
    const std::string dir_path = parent_dir + "/" + dir_name;
    if (directory_exists(dir_path)) {
        return true;
    }

    const std::string tar_path = parent_dir + "/" + dir_name + ".tar";
    const std::string zst_path = tar_path + ".zst";

    if (!file_exists(tar_path)) {
        if (!file_exists(zst_path)) {
            std::cerr << "Neither " << dir_path << " nor " << zst_path
                      << " exists; cannot restore " << dir_name << std::endl;
            return false;
        }
        std::ifstream in(zst_path, std::ios::binary | std::ios::ate);
        if (!in) {
            std::cerr << "cannot open " << zst_path << std::endl;
            return false;
        }
        uint64_t compressed_size = static_cast<uint64_t>(in.tellg());
        in.seekg(0, std::ios::beg);
        if (!zstd_decompress_stream_from_stream(in, compressed_size, tar_path)) {
            return false;
        }
    }

    std::string tar_cmd = "tar -C " + shell_quote(parent_dir) +
                          " -xf " + shell_quote(tar_path);
    return execute_system_command(tar_cmd, "tar extract failed: " + tar_path);
}

bool directory_exists_legacy_unused(const std::string& dir_path) {
    #ifdef _WIN32
        // 若需兼容Windows，可补充实现
        return false;
    #else
        struct stat st;
        return (stat(dir_path.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
    #endif
    }
 
 void reversed_complementary(std::string &seq_i, int read_len)
 {
	 int i = 0;
	 for(i = read_len - 1; i >= read_len/2; i--){
		 char tmp = seq_i[read_len - 1 - i];
		 seq_i[read_len - 1 - i] = charTchar[(int)seq_i[i] - 1];
		 seq_i[i] =  charTchar[(int)tmp - 1];
 
	 }
		 
	 //seq_i[read_len] = '\0';
 }
 
 void reversed_qual(std::string &seq_i, int read_len)
 {
	 int i = 0;
	 for(i = read_len - 1; i  >= read_len/2; i--){		
		 char tmp = seq_i[read_len - 1 - i];
		 seq_i[read_len - 1 - i] = seq_i[i];
		 seq_i[i] = tmp;
		 }
 }
 static inline void un_compact_sa_rst(uint64_t map_cpt, std::vector<BASIC_MAP_rst> &read_exact_mapping_rst, int break_read_length, int read_length){
	 int8_t is_rev = map_cpt & 0x1; map_cpt >>= 1;
	 uint32_t POS = map_cpt & POS_MASK; map_cpt >>= POS_OFFSET;
	 uint32_t wb_ID = map_cpt & 0x3ffffff;
	 is_rev = !is_rev;
	 
	 if(!is_rev){
		 if(POS >= WB_POS_MAX){
			 POS -= (break_read_length - 1);
			 //POS-=149;
			 if(POS > REF_POS_MAX || POS < WB_POS_MAX){
				 fprintf(stderr, "FATAL ERROR un_compact_sa_rst1, wb_ID %d @ POS %d\n ", wb_ID, POS);
				 return;  //check
			 }
		 }else{
			 if(POS < (break_read_length - 1)){
				 fprintf(stderr, "FATAL ERROR, un_compact_sa_rst2, wb_ID %d @ POS %d\n ", wb_ID, POS);
				 return; //check
			 }
			 POS -= (break_read_length-1);
			 //POS -= 149;
		 }
	 }else{
		 if(POS <(read_length-break_read_length)){
			 fprintf(stderr, "FATAL ERROR, un_compact_sa_rst3, wb_ID %d @ POS %d\n ", wb_ID, POS);
			 return;
		 }
		 POS -= (read_length-break_read_length);
		 //POS += (read_length -1);
	 }
	 
 
	 //fprintf(stderr, "\n emplace read_ID %d bg %ld, ed %ld, barcode %d",read_ID, read_bwt_128[read_end_ID].bwt_bg, read_bwt_128[read_end_ID].bwt_ed, barcode[read_ID]);
	 //fprintf(stderr, "is_rev %d, POS %d, wb_ID %d\n \n", is_rev, POS, wb_ID);
	 read_exact_mapping_rst.emplace_back();
	 read_exact_mapping_rst.back().set_exact_result(wb_ID, POS, is_rev, read_length);
 }
 
 
 void SAM_OUT::write_line(int s, uint flag, std::string  &cigar_p, std::string &md_p, kseq_t * read,
			 std::string &chr_id_p, uint32_t ref_p, int nm , int as,  std::string &mate_chr_id_p, uint32_t mate_ref_p, int rev, int read_length, int exc_){
	 exc = exc_;
	 read_p = read;
	 int sam_cross = ref_p + read_length - mate_ref_p;
	 std::string mate_chr = mate_chr_id_p;
 
	 std::string read_seq1 = std::string(read->seq.s);
	 std::string qual = std::string(read->qual.s);
	 
	 if(ref_p < mate_ref_p){
		 sam_cross = mate_ref_p - ref_p + read_length;
	 }
	 if(rev){
		 reversed_complementary(read_seq1, read_length);
		 reversed_qual(qual, read_length);
		 sam_cross = - sam_cross;
	 
	 }
 
	 //std::string rn = std::string(read->name.s);
	 // std::size_t pos = rn.find("/");      // position of "live" in str
	   // std::string read_name = rn.substr(0, pos); 
 
	 if (mate_chr_id_p == chr_id_p) mate_chr="=";
	 else sam_cross = 0;
	 
 
	 sam_line.clear();
	 sam_line = (std::string(std::string(read->name.s)+ "\t" + std::to_string(flag & 0Xff) + "\t"
				 +chr_id_p + "\t" + std::to_string(ref_p) 
				 + "\t" + std::to_string(s) + "\t" + cigar_p + "\t" + mate_chr + "\t" 
				 +  std::to_string(mate_ref_p) + "\t"
				 + std::to_string(sam_cross) + "\t" +read_seq1 + "\t"
				 + qual + "\t" + "NM:i:" + std::to_string(nm) + "\t"
				 + "MD:Z:" + md_p + "\t" + "AS:i:" + std::to_string(as) + "\t" + "RG:Z:readset1"
			 )
			 );
 }
 
 void test(const char* fn){
	 FILE* fp1 = fopen(fn, "rb");
	 Read_BWT *tmp = (Read_BWT *)xcalloc(2000,sizeof(Read_BWT));
	 fread(tmp, sizeof(Read_BWT), 2000, fp1);
	 for(int i= 0; i< 2000; i++){
		 fprintf(stderr, "%llu \n", tmp[i].bwt_k);
	 }
 
	 fclose(fp1);
	 
 }
 bool compare(const Compress_block& a, const Compress_block& b) {
	 if (a.window_id != b.window_id) {
		 return a.window_id < b.window_id;
	 }
	 return a.hap_id < b.hap_id;
 }
 bool create_directory(const std::string& dir_path) 
 {
	#ifdef _WIN32
		// Windows 系统：使用 _mkdir（返回 0 表示成功）
		if (_mkdir(dir_path.c_str()) == 0) {
			return true;
		}
		// 若文件夹已存在，也返回 true
		int err = errno;
		if (err == EEXIST) {
			return true;
		}
	#else
		// Linux/macOS 系统：使用 mkdir（权限 0755，返回 0 表示成功）
		if (mkdir(dir_path.c_str(), 0755) == 0) {
			return true;
		}
		// 若文件夹已存在，也返回 true
		int err = errno;
		if (err == EEXIST) {
			return true;
		}
	#endif
		// 创建失败
		std::cerr << "无法创建文件夹: " << dir_path << std::endl;
		return false;
 } 
std::string get_parent_dir(const char* pos_dir) {
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
 std::vector<unsigned char> convertStringToBits(const std::string& bit_str) {
    std::vector<unsigned char> buffer;
    size_t bit_count = bit_str.size();
    size_t byte_count = (bit_count + 7) / 8;  // 计算需要的字节数

    buffer.resize(byte_count, 0);  // 初始化缓冲区

    for (size_t i = 0; i < bit_count; ++i) {
        // 将字符'0'/'1'转换为位值（0/1）
        if (bit_str[i] == '1') {
            // 计算当前字节和位的位置
            size_t byte_idx = i / 8;
            int bit_offset = 7 - (i % 8);  // 高位优先存储
            buffer[byte_idx] |= (1 << bit_offset);
        }
    }
    return buffer;
}
 void store_the_quality_score(std::vector<std::string> quality_score_blocks, int id, char *quality_score_dir)
 {
	std::string filename = "quality_score." + std::to_string(id) + ".bin";
	filename = std::string(quality_score_dir) + '/' + filename; 
	std::cout<<"now write quality score to file:"<<filename<<std::endl;
	std::ofstream out_file(filename, std::ios::binary);
	if (!out_file) {
		 std::cerr << "无法打开文件: " << filename << std::endl;
		 exit(0);
	}
	size_t block_count = quality_score_blocks.size();
	out_file.write(reinterpret_cast<const char*>(&block_count), sizeof(size_t));
	for (std::string str : quality_score_blocks)
	{
		size_t code_len = str.length();
		out_file.write(reinterpret_cast<const char*>(&code_len), sizeof(size_t));
		out_file.write(str.data(), code_len);
	}
	std::cout<<"successfully write "<<filename<<std::endl;
 }
 bool store_the_final_block(const Compress_final_block_with_huffman_table& data, int id, char *pos_dir)
{
    std::string filename = "store." + std::to_string(id) + ".bin";
    filename = std::string(pos_dir) + '/' + filename; 
    // 以二进制模式打开文件
    std::cout << "now write final block to file:" << filename << std::endl;
    std::ofstream out_file(filename, std::ios::binary);
    if (!out_file) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        exit(0);
    }
    int block_count = data.Final_blocks.size();
    out_file.write(reinterpret_cast<const char*>(&block_count), sizeof(int));
    std::cout << block_count << std::endl;

    // 流水号（仅用于txt文件，从0开始）
    int serial_number = 0;

    for (const auto& block : data.Final_blocks) {
        // 原有二进制写入逻辑
        std::vector<unsigned char> binary_data = convertStringToBits(block.huffman_code);
        uint8_t code_len = block.huffman_code.size();
        int debug_code_len = block.huffman_code.size();
        if(debug_code_len>255)
        {
            std::cout<<"debug\n";
            std::cout<<"code_len: "<<code_len<<std::endl;
        }

        out_file.write(reinterpret_cast<const char*>(&code_len), sizeof(uint8_t));
        out_file.write(reinterpret_cast<const char*>(binary_data.data()), binary_data.size());
        /*
        // 写入qname元数据
        out_file.write(reinterpret_cast<const char*>(&block.lqname), sizeof(int));
        out_file.write(block.qname, block.lqname);
        
        // 写入序列1
        int real_seq1_len = block.real_seq1.size();
        out_file.write(reinterpret_cast<const char*>(&real_seq1_len), sizeof(int));
        out_file.write(block.real_seq1.data(), real_seq1_len);

        // 写入序列2
        int real_seq2_len = block.real_seq2.size();
        out_file.write(reinterpret_cast<const char*>(&real_seq2_len), sizeof(int));
        out_file.write(block.real_seq2.data(), real_seq2_len);
		*/

        //out_file.write(reinterpret_cast<const char*>(&block.hn), sizeof(uint8_t));
        //out_file.write(reinterpret_cast<const char*>(&block.hn2), sizeof(uint8_t));
    }

    // 写入哈夫曼表（原有逻辑）
    std::cout << "now write huffmantable" << std::endl;
    for (const auto& single_code : data.huffmanCode) {
        std::cout << single_code.first << ":" << single_code.second << std::endl;
        out_file.put(single_code.first);
        uint16_t code_len = single_code.second.size();
        out_file.write(reinterpret_cast<const char*>(&code_len), sizeof(uint16_t));
        out_file.write(single_code.second.data(), code_len);
    }

    std::cout << "successfully write " << filename << std::endl;
    return out_file.good();
}

void store_the_flag(const std::vector<std::pair<uint16_t, uint16_t>>& flag_blocks,
                    const char* pos_dir, int id)
{
    if (flag_blocks.empty()) {
        std::cerr << "flag_blocks为空，不写入文件" << std::endl;
        return;
    }

    std::string parent_dir = get_parent_dir(pos_dir);

    // ==================== 方案A：原始完整flag（/flags/） ====================
    {
        std::string flags_dir = parent_dir + "/flags";
        if (create_directory(flags_dir)) {
            std::string file_path = flags_dir + "/flags_" + std::to_string(id) + ".txt";
            std::ofstream out(file_path, std::ios::out | std::ios::trunc);
            if (!out.is_open()) {
                std::cerr << "无法打开文件(flags): " << file_path << std::endl;
            } else {
                int cnt = 0;
                for (const auto& pair : flag_blocks) {
                    out << cnt << " " << pair.first << " " << pair.second << "\n";
                    cnt++;
                }
                out.close();
                std::cout << "flags(完整)已写入: " << file_path
                          << "（记录数: " << flag_blocks.size() << "）" << std::endl;
            }
        }
    }

    // ==================== 方案B：仅保留反转位打包（/byte_flags/） ====================
    {
        std::string byte_flags_dir = parent_dir + "/byte_flags";
        if (create_directory(byte_flags_dir)) {
            std::string file_path = byte_flags_dir + "/flags_" + std::to_string(id) + ".bin";
            std::ofstream out(file_path, std::ios::binary | std::ios::trunc);
            if (!out.is_open()) {
                std::cerr << "无法打开文件(byte_flags): " << file_path << std::endl;
            } else {
                // 先写总对数
                uint64_t pair_count = flag_blocks.size();
                out.write(reinterpret_cast<const char*>(&pair_count), sizeof(pair_count));

                // 每对2bit：bit(2i)=R1反转, bit(2i+1)=R2反转，4对/字节
                const size_t byte_count = (flag_blocks.size() * 2 + 7) / 8;
                std::vector<uint8_t> packed(byte_count, 0);
                for (size_t i = 0; i < flag_blocks.size(); ++i) {
                    if (flag_blocks[i].first  & 0x10)
                        packed[(i * 2)     / 8] |= (1u << ((i * 2)     % 8));
                    if (flag_blocks[i].second & 0x10)
                        packed[(i * 2 + 1) / 8] |= (1u << ((i * 2 + 1) % 8));
                }

                out.write(reinterpret_cast<const char*>(packed.data()), packed.size());
                out.close();
                std::cout << "flags(压缩)已写入: " << file_path
                          << "（对数: " << pair_count
                          << "，字节数: " << byte_count << "）" << std::endl;
            }
        }
    }
}

 void store_the_name(const std::vector<std::string>& name_blocks, const char* pos_dir, int id) {
    if (name_blocks.empty()) {
        std::cerr << "name_blocks为空,不写入文件" << std::endl;
        return;
    }

    // 1. 构建目标文件路径：pos_dir的父路径 + "/flags.txt"
    std::string parent_dir = get_parent_dir(pos_dir);
	std::string flags_dir = parent_dir + "/qnames";
	if (!create_directory(flags_dir)) {
        return;  // 文件夹创建失败，直接返回
    }
    // 3. 构建文件名：flags/flags_id.txt（如 flags/flags_0.txt）
    std::string file_name = "qnames_" + std::to_string(id) + ".txt";
    std::string file_path = flags_dir + "/" + file_name;

    // 4. 打开文件并写入数据（竖排，每行一个pair）
    std::ofstream out(file_path, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "无法打开文件: " << file_path << std::endl;
        return;
    }
	int cnt = 0;
    // 5. 写入flag_blocks（每个pair占一行，两个元素用空格分隔）
    for (const auto& single : name_blocks) {
        out << cnt << " " << single << "\n";  // 竖排展示
		cnt++;
    }

    // 6. 关闭文件并输出提示
    out.close();
    std::cout << "name_blocks已写入: " << file_path 
              << "（记录数: " << name_blocks.size() << "）" << std::endl;
 }

 void store_the_diff_seq(const std::vector<std::string>& diff_seq_blocks, const char* pos_dir, int id) {
    if (diff_seq_blocks.empty()) {
        std::cout << "store_the_diff_seq: diff_seq_blocks为空，跳过写入。" << std::endl;
        return;
    }

    std::string parent_dir = get_parent_dir(pos_dir);
    std::string diff_seq_dir = parent_dir + "/diff_seq";
    if (!create_directory(diff_seq_dir)) {
        std::cerr << "store_the_diff_seq: 无法创建目录 " << diff_seq_dir << std::endl;
        return;
    }

    std::string file_path = diff_seq_dir + "/diff_" + std::to_string(id) + ".bin";
    std::ofstream out(file_path, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "store_the_diff_seq: 无法打开文件 " << file_path << std::endl;
        return;
    }

    const char magic[4] = {'D', 'S', 'E', 'Q'};
    uint64_t block_count = diff_seq_blocks.size();
    uint32_t read_length = static_cast<uint32_t>(g_global_read_length);
    uint64_t bit_count = block_count * read_length;
    uint64_t byte_count = (bit_count + 7) / 8;

    std::vector<uint8_t> packed(byte_count, 0);
    for (uint64_t row = 0; row < block_count; ++row) {
        const std::string& seq = diff_seq_blocks[row];
        if (seq.size() != read_length) {
            std::cerr << "store_the_diff_seq: 警告，ID " << id
                      << " 第 " << row << " 条长度为 " << seq.size()
                      << "，预期 " << read_length << std::endl;
        }
        const size_t usable = std::min<size_t>(seq.size(), read_length);
        for (size_t col = 0; col < usable; ++col) {
            if (seq[col] == '1') {
                uint64_t bit_idx = row * read_length + col;
                packed[bit_idx / 8] |= static_cast<uint8_t>(1u << (7 - (bit_idx % 8)));
            }
        }
    }

    out.write(magic, sizeof(magic));
    out.write(reinterpret_cast<const char*>(&block_count), sizeof(block_count));
    out.write(reinterpret_cast<const char*>(&read_length), sizeof(read_length));
    out.write(reinterpret_cast<const char*>(packed.data()), packed.size());
    out.close();

    if (out.good()) {
        std::cout << "store_the_diff_seq: 成功写入 " << file_path
                  << "（记录数: " << block_count
                  << "，字节数: " << byte_count << "）" << std::endl;
    } else {
        std::cerr << "store_the_diff_seq: 写入文件 " << file_path << " 时发生错误。" << std::endl;
    }
}

 void count_ACGT(std::vector<Compress_bio_string_block> bio_string_blocks, int id, char *pos_dir)
 {
	 Compress_final_block_with_huffman_table final_block_with_huffman_table;
	 std::unordered_map<char,uint64_t>freqMap;
	 for (const Compress_bio_string_block& block : bio_string_blocks)
	 {
		 int length = block.bio_string.size();
		 for (char c : block.bio_string) {
			 freqMap[c]++;
		 }
	 }
	 HuffmanCoder coder;
	 auto root = coder.buildHuffmanTree(freqMap);
	 coder.buildCodes(root, "", final_block_with_huffman_table.huffmanCode);
	 for (const Compress_bio_string_block& block : bio_string_blocks)
	 {
		 Compress_final_block final_block;
		 std::string encoded = coder.encode(block.bio_string, final_block_with_huffman_table.huffmanCode);
		 final_block.huffman_code = encoded;
		 int length_qname = strlen(block.qname);
		 memcpy(final_block.qname, block.qname, length_qname);
		 final_block.lqname = length_qname;
		 final_block.real_seq1 = block.real_seq1;
		 final_block.real_seq2 = block.real_seq2;
         final_block.hn = block.hn;
         final_block.hn2 = block.hn2;
         if(strcmp(block.qname,"E100003278L1C017R0114430563")==0)
         {
            std::cout<<"debug\n";
         }
		 final_block_with_huffman_table.Final_blocks.push_back(final_block);
	 }
	 if (store_the_final_block(final_block_with_huffman_table, id, pos_dir)) {
		 std::cout << "写入成功" << std::endl;
	 } else {
		 std::cerr << "写入失败" << std::endl;
	 }
 }

 void diff_seq_build(const std::vector<std::string>& diff_base_blocks, const char* pos_dir, int id) {
    if (diff_base_blocks.empty()) {
        std::cout << "diff_seq_build: 输入的 diff_base_blocks 为空，跳过处理。" << std::endl;
        return;
    }

    std::unordered_map<char, uint64_t> freqMap;
    for (const std::string& seq : diff_base_blocks) {
        for (char c : seq) {
            freqMap[c]++;
        }
        freqMap['\0']++; // 为每个序列添加结束符
    }
    if (freqMap.empty()) {
        std::cerr << "diff_seq_build: 没有统计到任何字符，无法构建哈夫曼树。" << std::endl;
        return;
    }

    // 2. 构建哈夫曼树和编码表
    HuffmanCoder coder;
    auto root = coder.buildHuffmanTree(freqMap);
    std::unordered_map<char, std::string> huffmanCode;
    coder.buildCodes(root, "", huffmanCode);

    // 3. 将所有碱基序列编码成一个完整的01字符串
    std::string full_binary_string;
    for (const std::string& seq : diff_base_blocks) {
        for (char c : seq) {
            full_binary_string += huffmanCode[c];
        }
        full_binary_string += huffmanCode['\0']; // 添加序列结束符的编码
    }

    if (full_binary_string.empty()) {
        std::cerr << "diff_seq_build: 编码后没有生成任何二进制数据，无法写入文件。" << std::endl;
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

    // 5. 准备文件路径并写入数据
    std::string parent_dir = get_parent_dir(pos_dir);
    std::string diff_seq_dir = parent_dir + "/diff_base";
    if (!create_directory(diff_seq_dir)) {
        std::cerr << "diff_seq_build: 无法创建目录 " << diff_seq_dir << std::endl;
        return;
    }

    std::string file_name = "diff_" + std::to_string(id) + ".bin";
    std::string file_path = diff_seq_dir + "/" + file_name;

    std::ofstream ofs(file_path, std::ios::binary);
    if (!ofs.is_open()) {
        std::cerr << "diff_seq_build: 无法打开文件 " << file_path << " 进行写入。" << std::endl;
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
        std::cout << "diff_seq_build: 成功将 " << (ofs.tellp()) << " 字节的数据写入到 " << file_path << std::endl;
    } else {
        std::cerr << "diff_seq_build: 写入文件 " << file_path << " 时发生错误。" << std::endl;
    }

    ofs.close();
}

 std::string num_to_ACGT(uint64_t num) {
    if (num == 0) return "A";
    std::string binary;
    uint64_t temp = num;
    while (temp > 0) {
        binary = (char)('0' + (temp & 1)) + binary;
        temp >>= 1;
    }
    if (binary.size() % 2 != 0) binary = "0" + binary;
    std::string result;
    for (size_t i = 0; i < binary.size(); i += 2) {
        char b0 = binary[i], b1 = binary[i+1];
        if      (b0=='0' && b1=='0') result += 'A';
        else if (b0=='0' && b1=='1') result += 'C';
        else if (b0=='1' && b1=='0') result += 'G';
        else                         result += 'T';
    }
    return result;
}

void sort_and_write_file(
        uint32_t &prev_window_id, int id,
        char *pos_dir, char *quality_score_dir, char *webp_dir,
        std::vector<Compress_block> blocks) {

    if (blocks.empty()) return;

    std::vector<std::string>                    quality_score_blocks;
    std::vector<Compress_bio_string_block>      bio_string_blocks;
    std::vector<std::pair<uint16_t, uint16_t>>  flag_blocks;
    std::vector<std::string>                    qname_blocks;
    std::vector<std::string>                    diff_seq_blocks;
    std::vector<std::string>                    diff_base_blocks;

    for (Compress_block& block : blocks) {
        const bool block_is_paired = block.is_paired != 0;
        
        quality_score_blocks.push_back(block.quality_score1);
        if (block_is_paired) quality_score_blocks.push_back(block.quality_score2);
        flag_blocks.push_back(std::make_pair(block.flags_1, block.flags_2));
        qname_blocks.push_back(std::string(block.qname));
        diff_seq_blocks.push_back(block.diff_seq1);
        if (block_is_paired) diff_seq_blocks.push_back(block.diff_seq2);
        diff_base_blocks.push_back(block.diff_base1);
        if (block_is_paired) diff_base_blocks.push_back(block.diff_base2);

        uint32_t window_id1 = (uint32_t)block.window_id;
        uint32_t window_id2 = (uint32_t)block.window_id2;

        // 字段1：window_dev（R1与上一R1的window_id差值）
        uint32_t window_dev = window_id1 - prev_window_id;

        // 字段3：pair_window_dev（同对内R2与R1的window_id差值，绝对值）
        uint32_t pair_window_dev = (window_id2 >= window_id1)
                                   ? (window_id2 - window_id1)
                                   : (window_id1 - window_id2);

        // 字段5：add_flag（R2 window 是否 ≥ R1）
        char add_flag = (window_id2 >= window_id1) ? 'A' : 'C';

        Compress_bio_string_block bio_string_block;
        // 字段1：window_dev
        bio_string_block.bio_string += num_to_ACGT((uint64_t)window_dev);
        bio_string_block.bio_string += '$';
        // 字段2：wb_ref_pos1
        bio_string_block.bio_string += num_to_ACGT((uint64_t)block.hap_offset);
        if (block_is_paired) {
            bio_string_block.bio_string += '$';
            // 字段3：pair_window_dev
            bio_string_block.bio_string += num_to_ACGT((uint64_t)pair_window_dev);
            bio_string_block.bio_string += '$';
            // 字段4：wb_ref_pos2
            bio_string_block.bio_string += num_to_ACGT((uint64_t)block.hap_offset2);
            bio_string_block.bio_string += '$';
            // 字段5：add_flag（最后字段，无$）
            bio_string_block.bio_string += add_flag;
        }
        memset(bio_string_block.qname, 0, sizeof(bio_string_block.qname));
        int length_qname = strlen(block.qname);
        memcpy(bio_string_block.qname, block.qname, length_qname);
        bio_string_block.lqname    = length_qname;
        bio_string_block.real_seq1 = block.real_seq1;
        bio_string_block.real_seq2 = block.real_seq2;
        bio_string_block.hn2 = block_is_paired ? 1 : 0;
        bio_string_blocks.push_back(std::move(bio_string_block));

        prev_window_id = window_id1; // 更新，供下一条使用
    }

    store_the_quality_score(quality_score_blocks, id, quality_score_dir);
    count_ACGT(bio_string_blocks, id, pos_dir);
    store_the_flag(flag_blocks, pos_dir, id);
    if (g_global_o && g_global_o->discard_qname) {
        std::cout << "[QNAME] discard_qname=1, skip qname stream for partition "
                  << id << std::endl;
    } else {
        store_the_name(qname_blocks, pos_dir, id);
    }
    store_the_diff_seq(diff_seq_blocks, pos_dir, id);
    diff_seq_build(diff_base_blocks, pos_dir, id);
    fprintf(stderr, " file size  %llu \n", (unsigned long long)blocks.size());
    blocks.clear();
}

bool read_compressed_block(Compress_final_block_with_huffman_table& blocks_with_huffman_table, std::ifstream& in_file) 
{
    // 读取Final_blocks数量
    int block_count;
    in_file.read(reinterpret_cast<char*>(&block_count), sizeof(int));
    blocks_with_huffman_table.Final_blocks.resize(block_count);
    std::cout << block_count << std::endl;
    // 读取每个压缩块
    for (auto& block : blocks_with_huffman_table.Final_blocks) 
    {
        // 读取binary_data的字节长度
        
        uint8_t code_len;
        if (!in_file.read(reinterpret_cast<char*>(&code_len), sizeof(uint8_t)) || 
            code_len > 1024 * 1024)  // 限制最大1MB防止溢出
        {
            std::cerr << "无效的binary_data长度" << std::endl;
            return false;
        }
        int binary_data_len = (code_len + 7)/8;
		std::vector<unsigned char> binary_data(binary_data_len);

        // 读取二进制数据
        in_file.read(reinterpret_cast<char*>(binary_data.data()), binary_data_len);
        if (!in_file) 
        {
            std::cerr << "读取binary_data失败" << std::endl;
            return false;
        }
        block.huffman_code = convertBitsToString(binary_data, code_len);
        /*
        if (!in_file.read(reinterpret_cast<char*>(&block.lqname), sizeof(int)) ||
            block.lqname < 0 || block.lqname > 1024 * 1024)  // 增加长度校验
        {
            std::cerr << "读取qname长度失败" << std::endl;
            return false;
        }

        // 读取qname数据
        in_file.read(&block.qname[0], block.lqname);
        if (!in_file) 
        {
            std::cerr << "读取qname内容失败" << std::endl;
            return false;
        }
        if(strcmp(block.qname,"E100003278L1C011R0371909053")==0)
         {
            std::cout<<"debug\n";
         }
		// 新增：读取real_seq1
		int real_seq1_len = 0;
		if (!in_file.read(reinterpret_cast<char*>(&real_seq1_len), sizeof(int))) {
    		std::cerr << "读取real_seq1长度失败" << std::endl;
    		return false;
		}
		block.real_seq1.resize(real_seq1_len);
		if (!in_file.read(&block.real_seq1[0], real_seq1_len)) {
    		std::cerr << "读取real_seq1内容失败" << std::endl;
    		return false;
		}

		// 新增：读取real_seq2
		int real_seq2_len = 0;
		if (!in_file.read(reinterpret_cast<char*>(&real_seq2_len), sizeof(int))) {
    		std::cerr << "读取real_seq2长度失败" << std::endl;
    		return false;
		}
		block.real_seq2.resize(real_seq2_len);
		if (!in_file.read(&block.real_seq2[0], real_seq2_len)) {
    		std::cerr << "读取real_seq2内容失败" << std::endl;
    		return false;
    	}
        */
	}

    // 读取完整的哈夫曼表（原代码固定读5个，改为读至文件结束）
    blocks_with_huffman_table.huffmanCode.clear();
    while (in_file.peek() != EOF)  // 改为动态读取所有表项
    {
        char ch;
        if (!in_file.get(ch)) {
            std::cerr << "读取字符失败" << std::endl;
            return false;
        }

        uint16_t code_len;
        if (!in_file.read(reinterpret_cast<char*>(&code_len), sizeof(uint16_t)) ||
            code_len > 1024 * 1024)  // 增加长度校验
        {
            std::cerr << "读取编码长度失败" << std::endl;
            return false;
        }

        std::string code(code_len, '\0');
        if (!in_file.read(&code[0], code_len)) {
            std::cerr << "读取编码内容失败" << std::endl;
            return false;
        }

        blocks_with_huffman_table.huffmanCode[ch] = code;
    }

    read_sum += blocks_with_huffman_table.Final_blocks.size();
    return true;
}

std::vector<std::string> webp_read(const char* pos_dir, int id) {
    std::vector<std::string> diff_seq_blocks;

    const int PIXELS_PER_SEQUENCE = g_global_read_length; // 必须与编码时的 pixels_per_sequence 保持一致

    std::string parent_dir = get_parent_dir(pos_dir);
    std::string webp_dir_path = parent_dir + "/diff_seq/combined_output." + std::to_string(id);

    std::vector<std::string> webp_files;
    DIR* dir = opendir(webp_dir_path.c_str());
    if (!dir) {
        std::cerr << "webp_read: 无法打开目录: " << webp_dir_path << std::endl;
        return diff_seq_blocks;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_type == DT_REG) {
            std::string filename(entry->d_name);
            if (filename.size() >= 5 && filename.compare(filename.size() - 5, 5, ".webp") == 0) {
                webp_files.push_back(webp_dir_path + "/" + filename);
            }
        }
    }
    closedir(dir);

    if (webp_files.empty()) {
        std::cerr << "webp_read: 在目录 " << webp_dir_path << " 中未找到任何 .webp 文件。" << std::endl;
        return diff_seq_blocks;
    }
    
    // 关键：必须对文件列表进行排序，以保证恢复顺序与写入顺序一致
    // 假设文件名是 part_0.webp, part_1.webp... 这种格式
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
            auto partial_matrix = RestoreMatrixFromWebpTiled(file_path, PIXELS_PER_SEQUENCE);
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

std::vector<std::string> read_diff_seq_blocks(const char* pos_dir, int id) {
    std::vector<std::string> diff_seq_blocks;
    std::string parent_dir = get_parent_dir(pos_dir);
    std::string file_path = parent_dir + "/diff_seq/diff_" + std::to_string(id) + ".bin";

    std::ifstream in(file_path, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "read_diff_seq_blocks: 无法打开文件 " << file_path << std::endl;
        return diff_seq_blocks;
    }

    char magic[4];
    if (!in.read(magic, sizeof(magic)) ||
        !(magic[0] == 'D' && magic[1] == 'S' && magic[2] == 'E' && magic[3] == 'Q')) {
        std::cerr << "read_diff_seq_blocks: 文件 " << file_path << " 不是有效的diff_seq文件。" << std::endl;
        return diff_seq_blocks;
    }

    uint64_t block_count = 0;
    uint32_t file_read_length = 0;
    if (!in.read(reinterpret_cast<char*>(&block_count), sizeof(block_count)) ||
        !in.read(reinterpret_cast<char*>(&file_read_length), sizeof(file_read_length))) {
        std::cerr << "read_diff_seq_blocks: 读取文件头失败: " << file_path << std::endl;
        return diff_seq_blocks;
    }

    if (file_read_length == 0) {
        std::cerr << "read_diff_seq_blocks: read_length为0: " << file_path << std::endl;
        return diff_seq_blocks;
    }
    if (g_global_read_length > 0 && file_read_length != static_cast<uint32_t>(g_global_read_length)) {
        std::cerr << "read_diff_seq_blocks: 警告，文件read_length=" << file_read_length
                  << "，当前参数read_length=" << g_global_read_length << std::endl;
    }

    uint64_t bit_count = block_count * file_read_length;
    uint64_t byte_count = (bit_count + 7) / 8;
    std::vector<uint8_t> packed(byte_count, 0);
    if (byte_count > 0 &&
        !in.read(reinterpret_cast<char*>(packed.data()), static_cast<std::streamsize>(packed.size()))) {
        std::cerr << "read_diff_seq_blocks: 读取压缩数据失败: " << file_path << std::endl;
        return {};
    }

    diff_seq_blocks.reserve(block_count);
    for (uint64_t row = 0; row < block_count; ++row) {
        std::string seq;
        seq.reserve(file_read_length);
        for (uint32_t col = 0; col < file_read_length; ++col) {
            uint64_t bit_idx = row * file_read_length + col;
            bool bit = (packed[bit_idx / 8] >> (7 - (bit_idx % 8))) & 1u;
            seq.push_back(bit ? '1' : '0');
        }
        diff_seq_blocks.push_back(std::move(seq));
    }

    std::cout << "read_diff_seq_blocks: 成功恢复 " << diff_seq_blocks.size()
              << " 个序列块。" << std::endl;
    return diff_seq_blocks;
}

std::string diff_seq_read(const char* pos_dir, int id) {

    // 1. 构建文件路径
    std::string parent_dir = get_parent_dir(pos_dir);
    std::string file_path = parent_dir + "/diff_base/diff_" + std::to_string(id) + ".bin";

    // 2. 打开文件
    std::ifstream ifs(file_path, std::ios::binary);
    if (!ifs.is_open()) {
        std::cerr << "diff_seq_read: 无法打开文件 " << file_path << " 进行读取。" << std::endl;
        return "";
    }

    // --- 读取文件头部 ---
    char magic[4];
    if (!ifs.read(magic, sizeof(magic)) || !(magic[0] == 'H' && magic[1] == 'u' && magic[2] == 'f' && magic[3] == 'f')) {
        std::cerr << "diff_seq_read: 文件 " << file_path << " 不是有效的哈夫曼压缩文件。" << std::endl;
        ifs.close();
        return "";
    }

    uint8_t code_table_size;
    if (!ifs.read(reinterpret_cast<char*>(&code_table_size), sizeof(code_table_size))) {
        std::cerr << "diff_seq_read: 读取编码表大小失败。" << std::endl;
        ifs.close();
        return "";
    }

    std::unordered_map<char, std::string> huffmanCodeTable;
    for (int i = 0; i < static_cast<int>(code_table_size); ++i) {
        char c; uint8_t code_len;
        if (!ifs.read(&c, sizeof(c)) || !ifs.read(reinterpret_cast<char*>(&code_len), sizeof(code_len))) {
            std::cerr << "diff_seq_read: 读取编码表条目失败。" << std::endl;
            ifs.close();
            return "";
        }
        std::string code(code_len, '\0');
        if (!ifs.read(&code[0], code_len)) {
            std::cerr << "diff_seq_read: 读取编码字符串失败。" << std::endl;
            ifs.close();
            return "";
        }
        huffmanCodeTable[c] = code;
    }

    uint64_t original_bits_length;
    if (!ifs.read(reinterpret_cast<char*>(&original_bits_length), sizeof(original_bits_length))) {
        std::cerr << "diff_seq_read: 读取原始比特长度失败。" << std::endl;
        ifs.close();
        return "";
    }
    // --- 头部读取完毕 ---

    // 3. 读取压缩数据字节流
    std::vector<uint8_t> byte_data((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    if (byte_data.empty() && original_bits_length > 0) {
        std::cerr << "diff_seq_read: 文件 " << file_path << " 中没有找到压缩数据。" << std::endl;
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
    HuffmanCoder coder;
    std::shared_ptr<HuffmanNode> root = coder.rebuildTreeFromCodeTable(huffmanCodeTable);
    if (!root) {
        std::cerr << "diff_seq_read: 无法从编码表重建哈夫曼树。" << std::endl;
        return "";
    }
    
    // !!! 核心修改：直接返回解码后的长字符串
    std::string decoded_string = coder.decode(binary_string, root);

    if (decoded_string.empty()) {
        std::cerr << "diff_seq_read: 解码失败，得到空字符串。" << std::endl;
    } else {
        std::cout << "diff_seq_read: 成功解码，得到 " << decoded_string.size() << " 个字符。" << std::endl;
    }

    return decoded_string;
}

 bool load_from_file(const std::string& filename, Compress_final_block_with_huffman_table& block) 
 {
	 std::ifstream in_file(filename, std::ios::binary);
	 if (!in_file) return false;
	 // 读取结构体数据
	 if (!read_compressed_block(block, in_file)) 
	 {
		 std::cerr << "读取数据失败" << std::endl;
		 return false;
	 }
	 return true;
 }
 std::string decode_huffman(const std::string& encoded_str,const std::unordered_map<std::string, char>& reverse_map) 
 {
	 std::string result;
	 std::string current_code;
 
	 for (char bit : encoded_str) {
		 current_code += bit;
		 
		 // 正确用法：检查迭代器有效性
		 const auto it = reverse_map.find(current_code);
		 if (it != reverse_map.end()) { 
			 result += it->second;
			 current_code.clear();
		 }
	 }
 
	 
	 if (!current_code.empty())
	 {
		 throw std::runtime_error("解码失败：存在未匹配的编码前缀 - " + current_code);
	 }
 
	 return result;
 }
 uint32_t quick_pow_32(int di, int zhi) {
	uint32_t sum = 1;
	while (zhi!= 0) {
		if (zhi & 1) {    //这里是zhi%2==1的意思
			sum *= di;
		}
		zhi >>= 1;    //这里是zhi/=2的意思，即向右移动一位在二进制中即除以2
		di = di * di;
	}
	return sum;
}
 uint64_t reverse_dev_num(std::string dev_string)
 {
	 uint64_t final = 0;
	 int length = dev_string.size();
	 std::string num_string;
	 for(int i = 0; i < length; i++)
	 {
		 if(dev_string[i] == 'A')num_string += "00";
		 else if(dev_string[i] == 'C')num_string += "01";
		 else if(dev_string[i] == 'G')num_string += "10";
		 else if(dev_string[i] == 'T')num_string += "11";
	 }
	 int num_length = num_string.size();
	 //std::cout<<"num string is"<<num_string<<std::endl;
	 for(int i = num_length - 1;i >= 0; i--)
	 {
		 if(num_string[i] == '1')final += quick_pow_32(2, num_length - i - 1);
	 }
	 return final;
 }
 uint64_t reverse_offset_num(std::string offset_string)
 {
	 //CAATA 0100001100 = 268
	 uint64_t final = 0;
	 int length = offset_string.size();
	 std::string num_string;
	 for(int i = 0; i < length; i++)
	 {
		 if(offset_string[i] == 'A')num_string += "00";
		 else if(offset_string[i] == 'C')num_string += "01";
		 else if(offset_string[i] == 'G')num_string += "10";
		 else if(offset_string[i] == 'T')num_string += "11";
	 }
	 //std::cout<<"num string is"<<num_string<<std::endl;
	 int num_length = num_string.size();
	 int base = 1;
	 for(int i = num_length - 1;i >= 0; i--)
	 {
		 if(num_string[i] == '1')final += base;
		 base *= 2;
	 }
	 return final;
 }

std::vector<Compress_block> bio_string_to_num(
    std::vector<Compress_bio_string_block> bio_string_blocks,
    std::vector<Compress_block> *origin_blocks,
    uint32_t pre_window_id)     // 上一条R1的window_id
{
    std::vector<Compress_block> blocks;

    for (const Compress_bio_string_block& block : bio_string_blocks) {
        // paired-end uses 5 fields; single-end uses 2 fields.
        std::vector<std::string> fields;
        std::string cur;
        
        for (char c : block.bio_string) {
            if (c == '$') { fields.push_back(cur); cur.clear(); }
            else          { cur += c; }
        }
        fields.push_back(cur); // 最后一个字段 add_flag（无$）

        if (fields.size() != 2 && fields.size() != 5) {
            fprintf(stderr, "bio_string_to_num: 字段数不是2或5，跳过\n");
            continue;
        }

        // ── 解码各字段 ──
        uint32_t window_dev      = (uint32_t)reverse_dev_num(fields[0]);
        uint32_t wb_ref_pos1     = (uint32_t)reverse_dev_num(fields[1]);
        uint32_t pair_window_dev = fields.size() == 5 ? (uint32_t)reverse_dev_num(fields[2]) : 0;
        uint32_t wb_ref_pos2     = fields.size() == 5 ? (uint32_t)reverse_dev_num(fields[3]) : 0;
        bool     add_flag        = (fields.size() == 5 && fields[4] == "A"); // true=window_id2>=window_id1

        // ── 重建 window_id ──
        uint32_t window_id1 = pre_window_id + window_dev;
        uint32_t window_id2 = fields.size() == 5
                              ? (add_flag
                              ? (window_id1 + pair_window_dev)
                              : (window_id1 - pair_window_dev))
                              : 0;

        // ── 填充 block ──
        Compress_block reverse_block;
        reverse_block.window_id   = window_id1;
        reverse_block.hap_offset = wb_ref_pos1;
        reverse_block.window_id2  = window_id2;
        reverse_block.hap_offset2 = wb_ref_pos2;
        reverse_block.is_paired = fields.size() == 5 ? 1 : 0;

        memcpy(reverse_block.qname, block.qname, block.lqname);
        reverse_block.lqname    = block.lqname;
        reverse_block.real_seq1 = block.real_seq1;
        reverse_block.real_seq2 = block.real_seq2;
        if(strcmp(block.qname,"A00900:594:HJVHNDSX7:1:1353:22444:9392")==0)
        {
            std::cout<<"debug\n";
        }
        pre_window_id = window_id1; // 更新，供下一条使用
        blocks.push_back(std::move(reverse_block));
    }

    return blocks;
}
int create_output_dir_c(char* webp_dir, char* output_dir, int output_dir_len) {
    if (webp_dir == NULL || output_dir == NULL || output_dir_len <= 0) {
        fprintf(stderr, "参数错误：webp_dir 或 output_dir 为空\n");
        return -1;
    }

    // 1. 提取 webp_dir 的父目录（dirname 会直接修改输入字符串，需先复制一份避免破坏原路径）
    char webp_dir_copy[256];  // 假设路径长度不超过 255，可根据实际场景调整
    strncpy(webp_dir_copy, webp_dir, sizeof(webp_dir_copy) - 1);
    webp_dir_copy[sizeof(webp_dir_copy) - 1] = '\0';  // 确保字符串结束符

    char* parent_dir = dirname(webp_dir_copy);  // 提取父目录（如 "/home/data/webp" → "/home/data"）
    if (parent_dir == NULL) {
        fprintf(stderr, "提取父目录失败\n");
        return -1;
    }

    // 2. 拼接 output_dir 完整路径（父目录 + "/output_dir"）
    int ret = snprintf(output_dir, output_dir_len, "%s/output_dir", parent_dir);
    if (ret >= output_dir_len) {
        fprintf(stderr, "output_dir 路径长度超过缓冲区限制\n");
        return -1;
    }

    // 3. 调用系统命令 mkdir -p 递归创建目录（-p 确保父目录不存在时也能创建，且目录已存在不报错）
    char mkdir_cmd[512];
    snprintf(mkdir_cmd, sizeof(mkdir_cmd), "mkdir -p \"%s\"", output_dir);
    int cmd_ret = system(mkdir_cmd);
    if (WEXITSTATUS(cmd_ret) != 0) {  // 检查命令执行结果
        fprintf(stderr, "创建 output_dir 失败：命令执行错误（路径：%s）\n", output_dir);
        return -1;
    }

    printf("output_dir 创建成功：%s\n", output_dir);
    return 0;
}

int delta=0;
std::string fix_the_index(char *hs, char *ts, char *seq, uint16_t *hap_offset, int i, uint8_t *hn, 
	uint8_t *tn, int *max_match_length, std::string hap_string_buff)
{
	std::string middle="";
	std::string reverse_seq = "";
	memset(hs,0,sizeof(hs));
	memset(ts,0,sizeof(ts));
	middle.clear();
    delta += (i - (int)(*hn)) + 1;            
    *hap_offset += (i - (int)(*hn)) + 1;  
    *max_match_length -= (i - (unsigned int)(*hn)) + 1; 
	*hn = i + 1; 
    
    memcpy(hs, seq, *hn);
    hs[*hn] = '\0';  
    memcpy(ts, seq + *hn + *max_match_length, *tn);
    ts[*tn] = '\0';  
    for (int k = *hap_offset; k < *hap_offset + *max_match_length; k++) {
        if (k >= 0 && k < (int)hap_string_buff.size()) {  // 避免越界
            middle += hap_string_buff[k];
        }
    }
    reverse_seq = std::string(hs) + middle + std::string(ts);
    return reverse_seq;
}

int compare_seq(char *seq, std::string reverse_seq)
{
	for(int i=0;i<strlen(seq);i++)
	{
		if(seq[i]=='N')continue;
		else if(seq[i]!=reverse_seq[i])
		{
			return 0;
		}
	}
	return 1;
}

uint32_t error_count = 0;
uint32_t fix_count = 0;
uint32_t hap_offset_error = 0;

int find_best_hap(std::vector<hap_t> hap_list, std::vector<Variant> var_list, char *seq,
    uint16_t *hap_id, std::string ref, uint pos, std::string qname, 
    int window_id, std::string &diff_seq, std::string &diff_base, int *best_begin_, int min_hamming_distance)
{
    Window_t c_window_item;
	c_window_item.hap_list = hap_list;
	c_window_item.var_list = var_list;
    std::string hap_string_buff;
    int hamming_distance  = min_hamming_distance;
    int best_hap = 0;
    int best_begin = 0;
	for(int i = 0; i < hap_list.size();i++)
    {
        hap_string_buff.clear();
		c_window_item.get_hap_string(i, ref, hap_string_buff);
        int j = -g_global_read_length;
        int hap_size = hap_string_buff.size();
        for(; j < hap_size; j++)
        {
            int hamming_test = 0;
            int k = 0;
            for(; k < g_global_read_length; k++)
            {
                if(j+k<0||j+k>=hap_size)
                {
                    hamming_test++;
                }
                else if(seq[k]!=hap_string_buff[j+k])hamming_test++;
                if(hamming_test > hamming_distance)break;
            }
            if(hamming_test <= hamming_distance)
            {
                best_hap = i;
                best_begin = j;
                hamming_distance = hamming_test;
            }
        }
    }
    *hap_id = best_hap;
    *best_begin_ = best_begin;
    hap_string_buff.clear();
	c_window_item.get_hap_string(best_hap, ref, hap_string_buff);
    int test_hamming = 0;
	for(int i=0;i<g_global_read_length;i++)
    {
        if(best_begin+i<0||best_begin+i>=hap_string_buff.size())
        {
            diff_seq += '1';
            diff_base += seq[i];
            test_hamming++;
        }
        else if(seq[i]==hap_string_buff[best_begin+i])diff_seq += '0';
        else
        {
            diff_seq += '1';
            diff_base += seq[i];
            test_hamming++;
        }
    }
    hamming_count[test_hamming]++;
    if(qname=="E100003278L1C019R0042180231")
    {
        std::cout<<"window_id: "<<window_id<<std::endl;
        std::cout<<"hap_id: "<<best_hap<<std::endl;
        std::cout<<"hap_offset: "<<best_begin<<std::endl;
        std::cout<<"hap_list.size: "<<hap_list.size()<<std::endl;
    }
	return hap_string_buff.size();
}

struct RawBAMData {
    bam1_t* bam_ptr;  // 原始BAM数据（消费者负责销毁）
    std::string work_dir;  // 工作目录
};

// 2. 任务包结构体（消费者解析后使用）
struct BAMTask {
    std::string read_name;
    int lqname;
    uint16_t flags;
    int32_t pos;
    uint local_offset;
    std::string seq;
    int chr_id;
    std::vector<uint32_t> cigar;
    std::string md;
    uint8_t md_len;
    std::vector<int> chr_bg_wb_ID;
    std::string quality_score;
    int error_flag;
};

// 2. 线程安全队列（生产者-消费者通信）
class SafeTaskQueue {
	private:
		std::queue<std::vector<RawBAMData>> queue_;
		std::mutex mutex_;
		std::condition_variable prod_cond_;
		std::condition_variable cons_cond_;
		std::atomic<bool> is_stop_;
		const size_t max_size_;
	
	public:
		SafeTaskQueue(size_t max_size) : is_stop_(false), max_size_(max_size) {}
	
		void push_batch(const std::vector<RawBAMData>& batch) {
			std::unique_lock<std::mutex> lock(mutex_);
			prod_cond_.wait(lock, [this]() { return queue_.size() < max_size_ || is_stop_; });
			if (is_stop_) return;
			queue_.push(batch);
			cons_cond_.notify_one();
		}
	
		bool pop_batch(std::vector<RawBAMData>& batch) {
			std::unique_lock<std::mutex> lock(mutex_);
			cons_cond_.wait(lock, [this]() { return !queue_.empty() || is_stop_; });
			if (is_stop_ && queue_.empty()) return false;
			batch = std::move(queue_.front());
			queue_.pop();
			prod_cond_.notify_one();
			return true;
		}
	
		void stop() {
			std::unique_lock<std::mutex> lock(mutex_);
			is_stop_ = true;
			prod_cond_.notify_all();
			cons_cond_.notify_all();
		}
	};

// 3. 线程安全结果集（存储最终Compress_block）
class SafeCompressBlocks {
	private:
		std::vector<Compress_block> blocks_;
		std::mutex mutex_;
		std::atomic<size_t> size_;  // 原子计数器：记录块数量（无需加锁即可读取）
	
	public:
		SafeCompressBlocks() : size_(0) {}  // 初始化大小为0
	
		// 添加块时，同时更新原子计数器
		void push(const Compress_block& block) {
			std::unique_lock<std::mutex> lock(mutex_);
			blocks_.push_back(block);
			size_.store(blocks_.size(), std::memory_order_relaxed);  // 非阻塞更新大小
		}
	
		// 获取当前大小（无锁，仅作粗略判断）
		size_t get_size() const {
			return size_.load(std::memory_order_relaxed);  // 无锁读取
		}
	
		// 加锁获取所有块（仅在需要处理时调用）
		std::vector<Compress_block> get_all() {
			std::unique_lock<std::mutex> lock(mutex_);
			return blocks_;
		}
	
		// 加锁清空并替换为新块
		void replace(const std::vector<Compress_block>& new_blocks) {
			std::unique_lock<std::mutex> lock(mutex_);
			blocks_.swap(const_cast<std::vector<Compress_block>&>(new_blocks));  // 高效替换
			size_.store(blocks_.size(), std::memory_order_relaxed);  // 更新大小
		}
	
		// 加锁清空所有块
		void clear() {
			std::unique_lock<std::mutex> lock(mutex_);
			blocks_.clear();
			size_.store(0, std::memory_order_relaxed);
		}
};

// 4. 全局共享变量（线程安全）
const size_t MAX_QUEUE_SIZE = 20;
SafeTaskQueue g_task_queue(MAX_QUEUE_SIZE);          // 任务队列（生产者→消费者）
std::mutex g_buffer_mutex;           // 保护buffer的互斥锁
std::atomic<uint64_t> g_success_reads{0};      // 成功进入主压缩 block 的 reads
std::atomic<uint64_t> g_error_reads{0};        // 写入 error.fastq 的 reads
std::atomic<uint64_t> g_producer_primary_reads{0};   // producer 读到的 primary reads
std::atomic<uint64_t> g_producer_paired_reads{0};    // producer 成功配对送入 queue 的 reads
std::atomic<uint64_t> g_producer_unpaired_reads{0};  // producer 没配上的 primary reads
std::atomic<uint64_t> g_skipped_secondary{0};
std::atomic<uint64_t> g_skipped_supplementary{0};
bool g_input_is_paired = true;
std::mutex g_findname_mutex;
CLASSIFY_SHARE_DATA* g_share = nullptr; 

Aln_online::IDX_loader* g_global_idx = nullptr; // 全局索引（只读，线程安全）
OL_PAR* g_global_o = nullptr;         // 全局参数（只读，线程安全）
const size_t BATCH_SIZE = 500000;
const size_t WRITE_THRESHOLD = 30000;
std::vector<std::unique_ptr<std::mutex>> g_partition_locks;
std::mutex g_fq_write_mutex; 
std::string g_error_fastq_path;
FILE* g_error_fq_fp = nullptr;
thread_local std::string tl_error_fq_buffer;
const size_t ERROR_FQ_FLUSH_THRESHOLD = 4 * 1024 * 1024;  // 4MB
static char g_error_fq_stdio_buffer[1 << 20];

void init_error_fastq_writer() {
    g_error_fastq_path = get_parent_dir(g_share->o->pos_dir) + "/error.fastq";
    g_error_fq_fp = fopen(g_error_fastq_path.c_str(), "w");
    if (!g_error_fq_fp) {
        fprintf(stderr, "Error: Failed to open error fq file %s\n", g_error_fastq_path.c_str());
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
    setvbuf(g_error_fq_fp, g_error_fq_stdio_buffer, _IOFBF, sizeof(g_error_fq_stdio_buffer));
}

void flush_error_fq_buffer_if_needed(bool force = false) {
    if (!force && tl_error_fq_buffer.size() < ERROR_FQ_FLUSH_THRESHOLD) {
        return;
    }
    if (tl_error_fq_buffer.empty()) {
        return;
    }

    std::lock_guard<std::mutex> lock(g_fq_write_mutex);
    if (!g_error_fq_fp) {
        fprintf(stderr, "Error: g_error_fq_fp is null when flushing error buffer\n");
        fflush(stderr);
        return;
    }

    size_t written = fwrite(tl_error_fq_buffer.data(), 1, tl_error_fq_buffer.size(), g_error_fq_fp);
    if (written != tl_error_fq_buffer.size()) {
        fprintf(stderr, "Error: Failed to fully write error fq buffer, expected=%zu actual=%zu\n",
                tl_error_fq_buffer.size(), written);
        fflush(stderr);
    }

    if (force) {
        fflush(g_error_fq_fp);
    }
    tl_error_fq_buffer.clear();
}

void close_error_fastq_writer() {
    if (g_error_fq_fp) {
        fflush(g_error_fq_fp);
        fclose(g_error_fq_fp);
        g_error_fq_fp = nullptr;
    }
}

void write_single_block_to_file(FILE* fp, const Compress_block& block) {
    // R1 固定成员
    fwrite(&block.window_id,   sizeof(block.window_id),   1, fp);
    fwrite(&block.hap_offset,  sizeof(block.hap_offset),  1, fp);  // wb_ref_pos1
    fwrite(block.qname,        sizeof(block.qname),        1, fp);
    fwrite(&block.lqname,      sizeof(block.lqname),       1, fp);
    fwrite(&block.flags_1,     sizeof(block.flags_1),      1, fp);

    // R1 字符串
    size_t q1_len = block.quality_score1.size();
    fwrite(&q1_len, sizeof(q1_len), 1, fp);
    fwrite(block.quality_score1.data(), q1_len, 1, fp);

    size_t s1_len = block.real_seq1.size();
    fwrite(&s1_len, sizeof(s1_len), 1, fp);
    fwrite(block.real_seq1.data(), s1_len, 1, fp);

    size_t d1_len = block.diff_seq1.size();
    fwrite(&d1_len, sizeof(d1_len), 1, fp);
    fwrite(block.diff_seq1.data(), d1_len, 1, fp);

    size_t db1_len = block.diff_base1.size();
    fwrite(&db1_len, sizeof(db1_len), 1, fp);
    fwrite(block.diff_base1.data(), db1_len, 1, fp);

    // R2 固定成员
    fwrite(&block.window_id2,  sizeof(block.window_id2),  1, fp);
    fwrite(&block.hap_offset2, sizeof(block.hap_offset2), 1, fp);  // wb_ref_pos2
    fwrite(block.qname2,       sizeof(block.qname2),       1, fp);
    fwrite(&block.lqname2,     sizeof(block.lqname2),      1, fp);
    fwrite(&block.flags_2,     sizeof(block.flags_2),      1, fp);
    fwrite(&block.is_paired,   sizeof(block.is_paired),    1, fp);

    // R2 字符串
    size_t q2_len = block.quality_score2.size();
    fwrite(&q2_len, sizeof(q2_len), 1, fp);
    fwrite(block.quality_score2.data(), q2_len, 1, fp);

    size_t s2_len = block.real_seq2.size();
    fwrite(&s2_len, sizeof(s2_len), 1, fp);
    fwrite(block.real_seq2.data(), s2_len, 1, fp);

    size_t d2_len = block.diff_seq2.size();
    fwrite(&d2_len, sizeof(d2_len), 1, fp);
    fwrite(block.diff_seq2.data(), d2_len, 1, fp);

    size_t db2_len = block.diff_base2.size();
    fwrite(&db2_len, sizeof(db2_len), 1, fp);
    fwrite(block.diff_base2.data(), db2_len, 1, fp);

    fwrite(&block.error_flag, sizeof(block.error_flag), 1, fp);
}

void flush_local_blocks(const std::vector<Compress_block>& blocks) {
    // 遍历所有block，按window_id分区写入
    for (const auto& block : blocks) {
        // 计算分区索引（window_id / 10000）
        int partition_idx = block.window_id / 10000;

        // 边界检查：确保分区索引在有效范围内
        if (partition_idx < 0 || partition_idx > chr_bg_wb_ID[genome.size()] / 10000) {
            continue;
        }

        // 获取对应分区的文件句柄和锁
        FILE* fp = g_share->blockFiles[partition_idx];
        if (!fp) {
            continue;
        }

        // 关键修改：加锁写入，确保同一分区同时只有一个线程写入
        std::lock_guard<std::mutex> lock(*g_partition_locks[partition_idx]);
        // 写入block（使用自定义序列化，避免直接写指针/string）
        write_single_block_to_file(fp, block);
    }
}

// 生产者线程：读取BAM文件，生成任务包
void producer_thread(htsFile* input_file, bam_hdr_t* header, char* work_dir) {
    bam1_t* read = bam_init1();
    std::vector<RawBAMData> current_batch;
    current_batch.reserve(BATCH_SIZE);

    bam1_t* r1 = nullptr;
    char* current_read_name = nullptr;

    ProducerPerf perf;
    bool is_first_in_batch = true;
    PerfClock::time_point batch_start;

    while (sam_read1(input_file, header, read) >= 0) 
    {
        if (read->core.flag & BAM_FSECONDARY) {
            g_skipped_secondary++;
            continue;
        }
        
        if (read->core.flag & BAM_FSUPPLEMENTARY) {
            g_skipped_supplementary++;
            continue;
        }
        
        g_producer_primary_reads++;
        char* read_name = bam_get_qname(read);

        if (is_first_in_batch) {
            batch_start = PerfClock::now();
            is_first_in_batch = false;
        }

        if (!g_input_is_paired) {
            bam1_t* single = bam_init1();
            bam_copy1(single, read);
            current_batch.push_back({single, std::string(work_dir)});
            g_producer_paired_reads++;
            perf.read_pairs++;

            if (current_batch.size() >= BATCH_SIZE) {
                auto before_push = PerfClock::now();
                perf.batch_fill_ms += elapsed_ms(batch_start, before_push);

                auto push_start = PerfClock::now();
                g_task_queue.push_batch(current_batch);
                auto push_end = PerfClock::now();
                perf.push_ms += elapsed_ms(push_start, push_end);

                perf.batch_cnt++;
                current_batch.clear();
                is_first_in_batch = true;
            }
            continue;
        }

        if (r1 == nullptr) {
            r1 = bam_init1();
            bam_copy1(r1, read);
            if (current_read_name) {
                free(current_read_name);
                current_read_name = nullptr;
            }
            current_read_name = strdup(read_name);
        } else {
            if (strcmp(read_name, current_read_name) == 0) {
                bam1_t* r2 = bam_init1();
                bam_copy1(r2, read);

                current_batch.push_back({r1, std::string(work_dir)});
                current_batch.push_back({r2, std::string(work_dir)});
                perf.read_pairs++;
                g_producer_paired_reads += 2;

                free(current_read_name);
                current_read_name = nullptr;
                r1 = nullptr;

                if (current_batch.size() >= BATCH_SIZE) {
                    auto before_push = PerfClock::now();
                    perf.batch_fill_ms += elapsed_ms(batch_start, before_push);

                    auto push_start = PerfClock::now();
                    g_task_queue.push_batch(current_batch);
                    auto push_end = PerfClock::now();
                    perf.push_ms += elapsed_ms(push_start, push_end);

                    perf.batch_cnt++;

                    if (perf.batch_cnt % 10 == 0) {
                        std::lock_guard<std::mutex> lock(g_perf_print_mutex);
                        fprintf(stderr,
                            "[producer] batches=%zu pairs=%zu avg_fill=%.2f ms avg_push=%.2f ms\n",
                            perf.batch_cnt,
                            perf.read_pairs,
                            perf.batch_fill_ms / perf.batch_cnt,
                            perf.push_ms / perf.batch_cnt);
                    }

                    current_batch.clear();
                    is_first_in_batch = true;
                }
            } else {
                g_producer_unpaired_reads++;
            
                fprintf(stderr,
                        "Warning: unpaired/non-adjacent primary read discarded in producer: %s\n",
                        current_read_name);
            
                bam_destroy1(r1);
                free(current_read_name);
                current_read_name = nullptr;
            
                r1 = bam_init1();
                bam_copy1(r1, read);
                current_read_name = strdup(read_name);
            }
        }
    }

    if (r1 != nullptr) {
        g_producer_unpaired_reads++;
    
        fprintf(stderr, "Warning: Unpaired read remains: %s\n", current_read_name);
        bam_destroy1(r1);
        r1 = nullptr;
        free(current_read_name);
        current_read_name = nullptr;
    }

    if (!current_batch.empty()) {
        auto before_push = PerfClock::now();
        perf.batch_fill_ms += elapsed_ms(batch_start, before_push);

        auto push_start = PerfClock::now();
        g_task_queue.push_batch(current_batch);
        auto push_end = PerfClock::now();
        perf.push_ms += elapsed_ms(push_start, push_end);

        perf.batch_cnt++;
    }

    g_task_queue.stop();
    bam_destroy1(read);

    if (perf.batch_cnt > 0) 
    {
        std::lock_guard<std::mutex> lock(g_perf_print_mutex);
        fprintf(stderr,
            "[producer-final] batches=%zu pairs=%zu avg_fill=%.2f ms avg_push=%.2f ms total_fill=%.2f s total_push=%.2f s\n",
            perf.batch_cnt,
            perf.read_pairs,
            perf.batch_fill_ms / perf.batch_cnt,
            perf.push_ms / perf.batch_cnt,
            perf.batch_fill_ms / 1000.0,
            perf.push_ms / 1000.0);
        fprintf(stderr,
            "[producer-final] batches=%zu pairs=%zu "
            "primary_reads=%llu paired_reads=%llu unpaired_reads=%llu "
            "secondary=%llu supplementary=%llu\n",
            perf.batch_cnt,
            perf.read_pairs,
            (unsigned long long)g_producer_primary_reads.load(),
            (unsigned long long)g_producer_paired_reads.load(),
            (unsigned long long)g_producer_unpaired_reads.load(),
            (unsigned long long)g_skipped_secondary.load(),
            (unsigned long long)g_skipped_supplementary.load());
        
    }
}

Compress_block return_error_block(bam1_t* read1, std::string qual_str1, bam1_t* read2, std::string qual_str2, 
        char* seq1, char* seq2)
{
    Compress_block block;
    const char* name1 = bam_get_qname(read1);
    const char* name2 = bam_get_qname(read2);
    block.window_id = chr_bg_wb_ID[genome.size()];
    block.window_id2 = chr_bg_wb_ID[genome.size()];

    memcpy(block.qname, name1, read1->core.l_qname);
    block.lqname = read1->core.l_qname;
    block.flags_1 = read1->core.flag;
    block.quality_score1 = qual_str1;
    block.real_seq1 = std::string(seq1);

    memcpy(block.qname2, name2, read2->core.l_qname);
    block.lqname2 = read2->core.l_qname;
    block.flags_2 = read2->core.flag;
    block.quality_score2 = qual_str2;
    block.real_seq2 = std::string(seq2);
    return block;
}

void write_error_reads_to_fq(const std::string& qname1, const std::string& seq1, const std::string& qual1,
                             const std::string& qname2, const std::string& seq2, const std::string& qual2) {
    if (tl_error_fq_buffer.capacity() == 0) {
        tl_error_fq_buffer.reserve(ERROR_FQ_FLUSH_THRESHOLD + 1024);
    }
    std::string out_qname1 = qname1;
    std::string out_qname2 = qname2;
    if (g_global_o && g_global_o->discard_qname) {
        uint64_t id = g_error_qname_counter.fetch_add(1);
        out_qname1 = "XZIP_ERROR_" + std::to_string(id);
        out_qname2 = out_qname1;
    }

    tl_error_fq_buffer += "@";
    tl_error_fq_buffer += out_qname1;
    tl_error_fq_buffer += " 1:N:0:\n";
    tl_error_fq_buffer += seq1;
    tl_error_fq_buffer += "\n+\n";
    tl_error_fq_buffer += qual1;
    tl_error_fq_buffer += "\n";

    tl_error_fq_buffer += "@";
    tl_error_fq_buffer += out_qname2;
    tl_error_fq_buffer += " 2:N:0:\n";
    tl_error_fq_buffer += seq2;
    tl_error_fq_buffer += "\n+\n";
    tl_error_fq_buffer += qual2;
    tl_error_fq_buffer += "\n";

    flush_error_fq_buffer_if_needed(false);
}

int get_cigar_head_clip_length(int cigar_n, const uint32_t *bam_cigar) {
    int head_clip_len = 0; // 开头剪切总长度（Soft+Hard Clip）

    for (int i = 0; i < cigar_n; ++i) {
        // 解析CIGAR长度和类型
        uint32_t length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
        int type = (int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK)); // 与原代码类型计算逻辑一致

        // 仅处理开头的Soft/Hard Clip，非剪切类型直接退出
        if (type == CIGAR_SOFT_CLIP || type == CIGAR_HARD_CLIP) {
            head_clip_len += length;
        } else {
            // 遇到第一个非剪切类型，终止循环（后续CIGAR与开头剪切无关）
            break;
        }
    }

    return head_clip_len;
}

bool build_read_block_fields(bam1_t* read,
                             char* seq,
                             std::string& qual_str,
                             int& window_id,
                             uint16_t& effective_offset,
                             std::string& diff_seq,
                             std::string& diff_base) {
    if (read->core.n_cigar == 0) {
        return false;
    }

    get_bam_seq(0, g_global_read_length, seq, read);
    uint8_t* qual = bam_get_qual(read);
    int qual_len = read->core.l_qseq;
    qual_str.clear();
    qual_str.reserve(qual_len);
    for (int j = 0; j < qual_len; ++j) {
        qual_str.push_back((char)(qual[j] + 33));
    }

    int chr_id = read->core.tid;
    if (chr_id < 0 || chr_id >= static_cast<int>(genome.size())) {
        return false;
    }

    uint wb_ref_pos = 0;
    std::string ref;
    const std::string& full_chrom_seq = genome[chr_id].second;
    uint local_offset = read->core.pos;
    if (local_offset == 0) window_id = chr_bg_wb_ID[chr_id];
    else window_id = chr_bg_wb_ID[chr_id] + (local_offset - 1) / g_global_read_length;
    local_offset -= (window_id - chr_bg_wb_ID[chr_id]) * g_global_read_length;
    if (local_offset >= g_global_read_length) {
        window_id++;
        local_offset -= g_global_read_length;
    }
    wb_ref_pos = local_offset;

    uint32_t st_pos = read->core.pos + 1 - wb_ref_pos;
    uint32_t start_idx = st_pos - 1;
    if (start_idx >= full_chrom_seq.size()) {
        return false;
    }
    uint32_t max_possible = full_chrom_seq.size() - start_idx;
    int actual_length = std::min(2 * g_global_read_length, (int)max_possible);
    if (actual_length <= 0) {
        return false;
    }
    ref.assign(full_chrom_seq, start_idx, actual_length);

    int head_clip = get_cigar_head_clip_length(read->core.n_cigar, bam_get_cigar(read));
    int eo = (int)wb_ref_pos - head_clip;
    if (eo < 0 && window_id > 0 && window_id > chr_bg_wb_ID[chr_id]) {
        window_id--;
        eo += g_global_read_length;

        uint32_t new_st_pos = (window_id - chr_bg_wb_ID[chr_id]) * g_global_read_length + 1;
        uint32_t new_start_idx = new_st_pos - 1;
        int new_len = std::min(2 * g_global_read_length,
                               (int)full_chrom_seq.size() - (int)new_start_idx);
        if (new_len <= 0) {
            return false;
        }
        ref.assign(full_chrom_seq, new_start_idx, new_len);
    } else if (eo < 0) {
        eo = 0;
    }
    effective_offset = (uint16_t)eo;

    diff_seq.clear();
    diff_base.clear();
    diff_seq.reserve(g_global_read_length);
    for (int k = 0; k < g_global_read_length; ++k) {
        int ref_k = eo + k;
        bool mismatch = (ref_k < 0 || ref_k >= (int)ref.size())
                        ? true
                        : (seq[k] != ref[ref_k]);
        if (mismatch) {
            diff_seq += '1';
            diff_base += seq[k];
        } else {
            diff_seq += '0';
        }
    }

    return true;
}

void write_error_read_to_fq(const std::string& qname,
                            const std::string& seq,
                            const std::string& qual,
                            int read_no) {
    if (tl_error_fq_buffer.capacity() == 0) {
        tl_error_fq_buffer.reserve(ERROR_FQ_FLUSH_THRESHOLD + 1024);
    }
    std::string out_qname = qname;
    if (g_global_o && g_global_o->discard_qname) {
        uint64_t id = g_error_qname_counter.fetch_add(1);
        out_qname = "XZIP_ERROR_" + std::to_string(id);
    }

    tl_error_fq_buffer += "@";
    tl_error_fq_buffer += out_qname;
    tl_error_fq_buffer += read_no == 2 ? " 2:N:0:\n" : " 1:N:0:\n";
    tl_error_fq_buffer += seq;
    tl_error_fq_buffer += "\n+\n";
    tl_error_fq_buffer += qual;
    tl_error_fq_buffer += "\n";

    flush_error_fq_buffer_if_needed(false);
}

// 7. 匹配相同read名称的逻辑
void consumer_thread() {
    thread_local std::vector<Compress_block> local_blocks;
    thread_local size_t block_count = 0;
    thread_local ConsumerPerf perf;

    local_blocks.reserve(WRITE_THRESHOLD);

    const size_t HEARTBEAT_EVERY_PAIRS = 10000;   // 每处理多少对 read 打一次心跳
    const size_t SUMMARY_EVERY_BATCHES = 1;       // 每多少个 batch 打一次汇总
    const double HEARTBEAT_EVERY_SEC = 5.0;       // 即使 pair 数没到，也每隔几秒报一次

    std::thread::id current_thread_id = std::this_thread::get_id();
    const size_t tid_hash = std::hash<std::thread::id>{}(current_thread_id);

    std::vector<RawBAMData> raw_batch;
    int total = 0;

    while (true) {
        auto pop_start = PerfClock::now();
        bool ok = g_task_queue.pop_batch(raw_batch);
        auto pop_end = PerfClock::now();
        perf.pop_wait_ms += elapsed_ms(pop_start, pop_end);

        if (!ok) {
            {
                std::lock_guard<std::mutex> lock(g_perf_print_mutex);
                fprintf(stderr,
                        "[consumer-exit %zu] queue stopped, exiting loop. thread_batches=%zu thread_pairs=%zu\n",
                        tid_hash, perf.batch_cnt, perf.pair_cnt);
                fflush(stderr);
            }
            break;
        }

        auto batch_start = PerfClock::now();

        if (raw_batch.empty()) {
            {
                std::lock_guard<std::mutex> lock(g_perf_print_mutex);
                fprintf(stderr, "[consumer-start %zu] got empty batch\n", tid_hash);
                fflush(stderr);
            }
            raw_batch.clear();
            continue;
        }

        if (g_input_is_paired && raw_batch.size() % 2 != 0) {
            fprintf(stderr, "Fatal error: Batch size is odd (not paired)\n");
            fflush(stderr);
            exit(EXIT_FAILURE);
        }

        const size_t reads_per_record = g_input_is_paired ? 2 : 1;
        const size_t total_pairs_in_batch = raw_batch.size() / reads_per_record;

        {
            std::lock_guard<std::mutex> lock(g_perf_print_mutex);
            fprintf(stderr,
                    "[consumer-start %zu] got batch: raw_batch.size()=%zu pairs=%zu, finished_thread_batches=%zu, pop_wait=%.2f ms\n",
                    tid_hash, raw_batch.size(), total_pairs_in_batch, perf.batch_cnt,
                    elapsed_ms(pop_start, pop_end));
            fflush(stderr);
        }

        auto last_heartbeat = PerfClock::now();

        for (size_t i = 0; i < raw_batch.size(); i += reads_per_record) {
            size_t pair_idx_in_batch = i / reads_per_record;

            // heartbeat：按 pair 数 or 按时间
            auto now_for_heartbeat = PerfClock::now();
            double sec_since_last_hb = std::chrono::duration<double>(now_for_heartbeat - last_heartbeat).count();
            if (pair_idx_in_batch == 0 ||
                (pair_idx_in_batch % HEARTBEAT_EVERY_PAIRS == 0) ||
                (sec_since_last_hb >= HEARTBEAT_EVERY_SEC)) {
                std::lock_guard<std::mutex> lock(g_perf_print_mutex);
                fprintf(stderr,
                    "[consumer-heartbeat %zu] batch_progress=%zu/%zu pairs, "
                    "local_blocks=%zu, block_count=%zu, "
                    "success_reads=%llu error_reads=%llu accounted_reads=%llu\n",
                    tid_hash,
                    pair_idx_in_batch,
                    total_pairs_in_batch,
                    local_blocks.size(),
                    block_count,
                    (unsigned long long)g_success_reads.load(),
                    (unsigned long long)g_error_reads.load(),
                    (unsigned long long)(g_success_reads.load() + g_error_reads.load()));
                fflush(stderr);
                last_heartbeat = now_for_heartbeat;
            }

            perf.pair_cnt++;

            bam1_t* read1 = raw_batch[i].bam_ptr;
            bam1_t* read2 = g_input_is_paired ? raw_batch[i + 1].bam_ptr : nullptr;

            const char* name1 = bam_get_qname(read1);
            const char* name2 = read2 ? bam_get_qname(read2) : "";

            Compress_block block;
            block.is_paired = g_input_is_paired ? 1 : 0;

            auto t_decode_0 = PerfClock::now();
            char seq1[400] = {0};
            char seq2[400] = {0};
            std::string qual_str1, qual_str2;
            int window_id1 = -1, window_id2 = -1;
            uint16_t effective_offset1 = 0, effective_offset2 = 0;
            std::string diff_seq1, diff_base1, diff_seq2, diff_base2;

            bool r1_valid = build_read_block_fields(read1, seq1, qual_str1, window_id1,
                                                    effective_offset1, diff_seq1, diff_base1);
            bool r2_valid = true;
            if (read2) {
                r2_valid = build_read_block_fields(read2, seq2, qual_str2, window_id2,
                                                   effective_offset2, diff_seq2, diff_base2);
            }
            auto t_decode_1 = PerfClock::now();
            perf.decode_ms += elapsed_ms(t_decode_0, t_decode_1);
            perf.ref_ms += elapsed_ms(t_decode_0, t_decode_1);
            perf.hap_ms += elapsed_ms(t_decode_0, t_decode_1);

            if (!r1_valid || !r2_valid) {
                auto t_err_0 = PerfClock::now();
                if (read2) {
                    write_error_reads_to_fq(name1, seq1, qual_str1, name2, seq2, qual_str2);
                    g_error_reads += 2;
                    total += 2;
                } else {
                    write_error_read_to_fq(name1, seq1, qual_str1, 1);
                    g_error_reads += 1;
                    total += 1;
                }
                auto t_err_1 = PerfClock::now();
                perf.error_io_ms += elapsed_ms(t_err_0, t_err_1);
                perf.error_pair_cnt++;

                auto t_clean_0 = PerfClock::now();
                bam_destroy1(read1);
                if (read2) bam_destroy1(read2);
                auto t_clean_1 = PerfClock::now();
                perf.cleanup_ms += elapsed_ms(t_clean_0, t_clean_1);
                continue;
            }

            auto t_block_0 = PerfClock::now();

            block.window_id = window_id1;
            block.hap_offset = effective_offset1;
            memcpy(block.qname, name1, read1->core.l_qname);
            block.lqname = read1->core.l_qname;
            block.flags_1 = read1->core.flag;
            block.quality_score1 = std::move(qual_str1);
            block.real_seq1 = std::string(seq1);
            block.diff_seq1 = std::move(diff_seq1);
            block.diff_base1 = std::move(diff_base1);

            if (read2) {
                block.window_id2 = window_id2;
                block.hap_offset2 = effective_offset2;
                memcpy(block.qname2, name2, read2->core.l_qname);
                block.lqname2 = read2->core.l_qname;
                block.flags_2 = read2->core.flag;
                block.quality_score2 = std::move(qual_str2);
                block.real_seq2 = std::string(seq2);
                block.diff_seq2 = std::move(diff_seq2);
                block.diff_base2 = std::move(diff_base2);
            }
            local_blocks.push_back(std::move(block));

            auto t_block_1 = PerfClock::now();
            perf.block_ms += elapsed_ms(t_block_0, t_block_1);

            block_count++;
            total += reads_per_record;
            g_success_reads += reads_per_record;

            if (block_count >= WRITE_THRESHOLD) {
                {
                    std::lock_guard<std::mutex> lock(g_perf_print_mutex);
                    fprintf(stderr,
                            "[consumer-flush-start %zu] block_count=%zu local_blocks=%zu\n",
                            tid_hash, block_count, local_blocks.size());
                    fflush(stderr);
                }

                auto t_flush_0 = PerfClock::now();
                flush_local_blocks(local_blocks);
                auto t_flush_1 = PerfClock::now();
                perf.flush_ms += elapsed_ms(t_flush_0, t_flush_1);
                perf.flush_cnt++;

                {
                    std::lock_guard<std::mutex> lock(g_perf_print_mutex);
                    fprintf(stderr,
                            "[consumer-flush-end %zu] flush_time=%.2f ms flush_cnt=%zu\n",
                            tid_hash, elapsed_ms(t_flush_0, t_flush_1), perf.flush_cnt);
                    fflush(stderr);
                }

                local_blocks.clear();
                block_count = 0;
            }

            auto t_clean_0 = PerfClock::now();
            bam_destroy1(read1);
            if (read2) bam_destroy1(read2);
            auto t_clean_1 = PerfClock::now();
            perf.cleanup_ms += elapsed_ms(t_clean_0, t_clean_1);
        }

        auto batch_end = PerfClock::now();
        perf.total_batch_ms += elapsed_ms(batch_start, batch_end);
        perf.batch_cnt++;

        if (perf.batch_cnt % SUMMARY_EVERY_BATCHES == 0) {
            std::lock_guard<std::mutex> lock(g_perf_print_mutex);
            double b = (double)perf.batch_cnt;
            double p = (double)std::max<size_t>(perf.pair_cnt, 1);

            fprintf(stderr,
                    "[consumer-summary %zu] batches=%zu pairs=%zu "
                    "avg_pop=%.2f ms/batch avg_total=%.2f ms/batch "
                    "decode=%.4f ms/pair ref=%.4f ms/pair error_io=%.4f ms/pair "
                    "hamming=%.4f ms/pair hap=%.4f ms/pair block=%.4f ms/pair "
                    "flush=%.4f ms/pair cleanup=%.4f ms/pair error_pairs=%zu flush_cnt=%zu total_reads=%d\n",
                    tid_hash,
                    perf.batch_cnt, perf.pair_cnt,
                    perf.pop_wait_ms / b,
                    perf.total_batch_ms / b,
                    perf.decode_ms / p,
                    perf.ref_ms / p,
                    perf.error_io_ms / p,
                    perf.hamming_ms / p,
                    perf.hap_ms / p,
                    perf.block_ms / p,
                    perf.flush_ms / p,
                    perf.cleanup_ms / p,
                    perf.error_pair_cnt,
                    perf.flush_cnt,
                    (unsigned long long)(g_success_reads.load() + g_error_reads.load()));
            fflush(stderr);
        }

        raw_batch.clear();
    }

    flush_error_fq_buffer_if_needed(true);

    if (!local_blocks.empty()) {
        {
            std::lock_guard<std::mutex> lock(g_perf_print_mutex);
            fprintf(stderr,
                    "[consumer-final-flush-start %zu] remaining_local_blocks=%zu\n",
                    tid_hash, local_blocks.size());
            fflush(stderr);
        }

        auto t_flush_0 = PerfClock::now();
        flush_local_blocks(local_blocks);
        auto t_flush_1 = PerfClock::now();
        perf.flush_ms += elapsed_ms(t_flush_0, t_flush_1);
        perf.flush_cnt++;

        {
            std::lock_guard<std::mutex> lock(g_perf_print_mutex);
            fprintf(stderr,
                    "[consumer-final-flush-end %zu] flush_time=%.2f ms flush_cnt=%zu\n",
                    tid_hash, elapsed_ms(t_flush_0, t_flush_1), perf.flush_cnt);
            fflush(stderr);
        }

        local_blocks.clear();
    }

    if (perf.batch_cnt > 0) {
        std::lock_guard<std::mutex> lock(g_perf_print_mutex);
        double b = (double)perf.batch_cnt;
        double p = (double)std::max<size_t>(perf.pair_cnt, 1);

        fprintf(stderr,
                "[consumer-final %zu] batches=%zu pairs=%zu "
                "avg_pop=%.2f ms/batch avg_total=%.2f ms/batch "
                "decode=%.4f ms/pair ref=%.4f ms/pair error_io=%.4f ms/pair "
                "hamming=%.4f ms/pair hap=%.4f ms/pair block=%.4f ms/pair "
                "flush=%.4f ms/pair cleanup=%.4f ms/pair error_pairs=%zu flush_cnt=%zu total_reads=%d\n",
                tid_hash,
                perf.batch_cnt, perf.pair_cnt,
                perf.pop_wait_ms / b,
                perf.total_batch_ms / b,
                perf.decode_ms / p,
                perf.ref_ms / p,
                perf.error_io_ms / p,
                perf.hamming_ms / p,
                perf.hap_ms / p,
                perf.block_ms / p,
                perf.flush_ms / p,
                perf.cleanup_ms / p,
                perf.error_pair_cnt,
                perf.flush_cnt,
                (unsigned long long)(g_success_reads.load() + g_error_reads.load()));
        fflush(stderr);
    }
}

bool compare_blocks(const Compress_block& a, const Compress_block& b) {
    if (a.window_id != b.window_id)
        return a.window_id < b.window_id;
    return a.hap_offset < b.hap_offset;  // 同 window 内按 wb_ref_pos 排序
}

// 读取单个分区文件中的所有Compress_block（严格匹配write_single_block_to_file的格式）
std::vector<Compress_block> read_compress_blocks(const std::string& file_path) {
    std::vector<Compress_block> blocks;
    FILE* fp = fopen(file_path.c_str(), "rb");
    if (!fp) {
        fprintf(stderr, "Warning: 无法打开分区文件 %s\n", file_path.c_str());
        return blocks;
    }

    // 循环读取所有block，直到文件结束
    while (true) {
        Compress_block block;
        size_t total_read = 0;  // 统计当前block的读取字节数，用于判断是否读到文件末尾
        // ==============================================
        // 读取R1固定大小成员（与写入顺序完全一致）
        // ==============================================
        total_read += fread(&block.window_id,   sizeof(block.window_id),   1, fp);
        total_read += fread(&block.hap_offset,  sizeof(block.hap_offset),  1, fp);  // wb_ref_pos1
        total_read += fread(block.qname,        sizeof(block.qname),        1, fp);
        total_read += fread(&block.lqname,      sizeof(block.lqname),       1, fp);
        total_read += fread(&block.flags_1,     sizeof(block.flags_1),      1, fp);
        

        // ==============================================
        // 读取R1的std::string成员（先长度后内容）
        // ==============================================
        size_t q1_len = 0;
        total_read += fread(&q1_len, sizeof(q1_len), 1, fp);
        block.quality_score1.resize(q1_len);
        if (q1_len > 0) {
            total_read += fread(block.quality_score1.data(), q1_len, 1, fp);
        }

        size_t s1_len = 0;
        total_read += fread(&s1_len, sizeof(s1_len), 1, fp);
        block.real_seq1.resize(s1_len);
        if (s1_len > 0) {
            total_read += fread(block.real_seq1.data(), s1_len, 1, fp);
        }

        size_t d1_len = 0;
        total_read += fread(&d1_len, sizeof(d1_len), 1, fp);
        block.diff_seq1.resize(d1_len);
        if (d1_len > 0) {
            total_read += fread(block.diff_seq1.data(), d1_len, 1, fp);
        }

        size_t db1_len = 0;
        total_read += fread(&db1_len, sizeof(db1_len), 1, fp);
        block.diff_base1.resize(db1_len);
        if (db1_len > 0) {
            total_read += fread(block.diff_base1.data(), db1_len, 1, fp);
        }

        // ==============================================
        // 读取R2固定大小成员（与写入顺序完全一致）
        // ==============================================
        total_read += fread(&block.window_id2,  sizeof(block.window_id2),  1, fp);
        total_read += fread(&block.hap_offset2, sizeof(block.hap_offset2), 1, fp);  // wb_ref_pos2
        total_read += fread(block.qname2,       sizeof(block.qname2),       1, fp);
        total_read += fread(&block.lqname2,     sizeof(block.lqname2),      1, fp);
        total_read += fread(&block.flags_2,     sizeof(block.flags_2),      1, fp);
        total_read += fread(&block.is_paired,   sizeof(block.is_paired),    1, fp);
        

        // ==============================================
        // 读取R2的std::string成员（先长度后内容）
        // ==============================================
        size_t q2_len = 0;
        total_read += fread(&q2_len, sizeof(q2_len), 1, fp);
        block.quality_score2.resize(q2_len);
        if (q2_len > 0) {
            total_read += fread(block.quality_score2.data(), q2_len, 1, fp);
        }

        size_t s2_len = 0;
        total_read += fread(&s2_len, sizeof(s2_len), 1, fp);
        block.real_seq2.resize(s2_len);
        if (s2_len > 0) {
            total_read += fread(block.real_seq2.data(), s2_len, 1, fp);
        }

        size_t d2_len = 0;
        total_read += fread(&d2_len, sizeof(d2_len), 1, fp);
        block.diff_seq2.resize(d2_len);
        if (d2_len > 0) {
            total_read += fread(block.diff_seq2.data(), d2_len, 1, fp);
        }

        size_t db2_len = 0;
        total_read += fread(&db2_len, sizeof(db2_len), 1, fp);
        block.diff_base2.resize(db2_len);
        if (db2_len > 0) {
            total_read += fread(block.diff_base2.data(), db2_len, 1, fp);
        }
        total_read += fread(&block.error_flag, sizeof(block.error_flag), 1, fp);
        if (total_read == 0) {
            break;
        }
        if(strcmp(block.qname,"E100003278L1C011R0371909053")==0)
         {
            std::cout<<"debug\n";
         }
        blocks.push_back(block);
    }

    fclose(fp);
    return blocks;
}

// 将排序后的blocks写回文件（严格复用你的write_single_block_to_file函数）
void write_compress_blocks(const std::string& file_path, const std::vector<Compress_block>& blocks) {
    FILE* fp = fopen(file_path.c_str(), "wb");  // "wb"模式覆盖原文件
    if (!fp) {
        fprintf(stderr, "Error: 无法写入分区文件 %s\n", file_path.c_str());
        return;
    }

    // 逐个写入block，直接复用你的写入函数，确保格式一致
    for (const auto& block : blocks) {
        write_single_block_to_file(fp, block);
    }

    fclose(fp);
}

void parallel_sort_and_write(
    std::vector<uint32_t>& pre_window_id_array,
    int total_groups,
    char* pos_dir, 
    char* quality_dir,
    char* webp_dir,
    const std::string& work_dir,
    int max_threads
) {
    std::vector<std::thread> threads;
    int active_threads = 0;

    for (int i = 0; i < total_groups; ++i) 
    {
        if (active_threads >= max_threads) {
            for (auto it = threads.begin(); it != threads.end(); ++it) {
                if (it->joinable()) {
                    it->join();
                    threads.erase(it);
                    active_threads--;
                    break;
                }
            }
        }

        threads.emplace_back([=]() {
            const std::string file_name = work_dir + "/tmp/test." + std::to_string(i) + ".bin";

            std::vector<Compress_block> blocks = read_compress_blocks(file_name);
            if (blocks.empty()) {
                fprintf(stderr, "Warning: No blocks in partition %d, skipping\n", i);
                return;
            }

            uint32_t prev_window_id = pre_window_id_array[i];
            sort_and_write_file(prev_window_id, i, pos_dir, quality_dir, webp_dir, blocks);
        });
        active_threads++;
    }

    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
}

void process_window_range(
    int chr_id,
    int start_window,
    int end_window,
    int window_size,
    int read_length,
    const std::string& full_chrom_seq,
    const std::vector<int>& chr_bg_wb_ID
) {
    try {
        for (int window_id = start_window; window_id <= end_window; window_id++) {
            Window_t& c_window = g_window_info[window_id];
            c_window.chr_ID = chr_id;
            c_window.st_pos = (window_id - chr_bg_wb_ID[chr_id]) * window_size + 1;

            uint32_t st_pos = c_window.st_pos;
            int length = read_length;
            uint32_t start_idx = st_pos - 1;
            uint32_t max_possible = full_chrom_seq.size() - start_idx;
            int actual_length = std::min(length, (int)max_possible);
            std::string ref;
            if (actual_length > 0) {
                ref.assign(full_chrom_seq, start_idx, actual_length);
            }
            std::vector<hap_t> hap_list(1);
            std::vector<Variant> var_list;
            hap_list[0].hap_length = ref.size();
            c_window.hap_list = hap_list;
            c_window.var_list = var_list;
        }
    } catch (...) {
        std::cerr << "线程 " << std::this_thread::get_id() << " 处理染色体 chr_id=" << chr_id << " 时出错" << std::endl;
    }

    std::cout << "线程 " << std::this_thread::get_id() << " 完成染色体 chr_id=" << chr_id << " 处理（窗口范围: " << start_window << "~" << end_window << "）" << std::endl;
}

void process_partition(
    int partition_idx,
    const std::string& work_dir,
    std::vector<uint32_t>& last_window_id,  // 改名，去掉 hap_id_list
    std::mutex& mtx
) {
    const std::string file_name = work_dir + "/tmp/test." + std::to_string(partition_idx) + ".bin";
    std::vector<Compress_block> blocks = read_compress_blocks(file_name);
    if (blocks.empty()) return;
    std::sort(blocks.begin(), blocks.end(), compare_blocks);
    last_window_id[partition_idx] = blocks.back().window_id;
    write_compress_blocks(file_name, blocks);
    fprintf(stderr, "Thread %d: Written sorted blocks to %s\n", partition_idx, file_name.c_str());
}

FastaData read_fasta(const std::string& fasta_path, int start_chr = -1, int end_chr = -1) {
    FastaData fasta_data;
    std::ifstream in(fasta_path);

    if (!in.is_open()) {
        std::cerr << "Error: 无法打开FASTA文件 " << fasta_path << std::endl;
        return fasta_data;
    }

    std::string line;
    std::string current_id;
    std::string current_seq;

    int chr_count = 0;
    bool keep_current = false;

    auto save_current = [&]() {
        if (!current_id.empty() && keep_current) {
            fasta_data.emplace_back(current_id, current_seq);

            std::cout << "ID: " << current_id
                      << ", 长度: " << current_seq.size()
                      << ", 前100个碱基: " << current_seq.substr(0, 100)
                      << std::endl;
        }

        current_id.clear();
        current_seq.clear();
        keep_current = false;
    };

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // 保存上一条染色体
            save_current();

            // 新染色体编号
            chr_count++;

            // 如果已经超过范围，直接退出
            if (end_chr > 0 && chr_count > end_chr) {
                break;
            }

            // 判断当前染色体是否需要读取
            keep_current = true;

            if (start_chr > 0 && chr_count < start_chr) {
                keep_current = false;
            }

            size_t space_pos = line.find(' ');
            current_id = (space_pos != std::string::npos)
                       ? line.substr(1, space_pos - 1)
                       : line.substr(1);

            current_seq.clear();
        } else {
            if (keep_current) {
                std::string upper_line = line;

                std::transform(
                    upper_line.begin(),
                    upper_line.end(),
                    upper_line.begin(),
                    [](unsigned char c) {
                        return std::toupper(c);
                    }
                );

                current_seq += upper_line;
            }
        }
    }

    // 保存最后一条
    save_current();

    in.close();
    return fasta_data;
}

std::vector<Window_t> generate_window_info(std::vector<int> chr_bg_wb_ID, int thread_num) 
{
    int total_chrs = genome.size();  // 总染色体数
    std::vector<std::pair<int, int>> chr_window_ranges;  // 每个染色体的窗口范围（start, end）

    // 1. 预处理每个染色体的窗口范围（不变）
    for (int chr_id = 0; chr_id < total_chrs; chr_id++) {
        int start = chr_bg_wb_ID[chr_id];
        int end = chr_bg_wb_ID[chr_id + 1] - 1;
        chr_window_ranges.emplace_back(start, end);
    }

    // 2. 预分配全局窗口数组（不变）
    int max_window_id = 0;
    for (const auto& range : chr_window_ranges) {
        max_window_id = std::max(max_window_id, range.second);
    }
    g_window_info.resize(max_window_id + 1);

    // 3. 多线程处理：一个线程负责一整个染色体（核心修改）
    std::vector<std::thread> threads;
    int window_size = g_global_read_length;
    int read_length = 2 * g_global_read_length;

    // 优化：线程数不超过染色体数（避免创建无用线程）
    int actual_thread_num = std::min(thread_num, total_chrs);
    std::cout << "实际使用线程数: " << actual_thread_num << "（总染色体数: " << total_chrs << "）" << std::endl;

    // 3.1 先启动第一批线程（填满线程池）
    int chr_id = 0;
    for (; chr_id < actual_thread_num; chr_id++) {
        const auto& [start_window, end_window] = chr_window_ranges[chr_id];
        if (start_window > end_window) {
            std::cout << "跳过空染色体 chr_id: " << chr_id << std::endl;
            continue;
        }

        const auto& chr_pair = genome[chr_id];
        const std::string& full_chrom_seq = chr_pair.second;

        std::cout << "线程处理染色体 chr_id: " << chr_id << "（窗口范围: " << start_window << "~" << end_window << "）" << std::endl;

        // 启动线程：传递整个染色体的窗口范围（start_window ~ end_window）
        threads.emplace_back(
            process_window_range,
            chr_id,
            start_window,  // 整个染色体的起始窗口
            end_window,    // 整个染色体的结束窗口（不拆分）
            window_size,
            read_length,
            full_chrom_seq,
            chr_bg_wb_ID
        );
    }

    // 3.2 等待第一批线程完成后，处理剩余染色体（避免线程数过多）
    for (; chr_id < total_chrs; chr_id++) {
        // 等待一个线程完成，释放资源
        if (!threads.empty()) {
            threads.back().join();
            threads.pop_back();
        }

        const auto& [start_window, end_window] = chr_window_ranges[chr_id];
        if (start_window > end_window) {
            std::cout << "跳过空染色体 chr_id: " << chr_id << std::endl;
            continue;
        }

        const auto& chr_pair = genome[chr_id];
        const std::string& full_chrom_seq = chr_pair.second;

        std::cout << "线程处理染色体 chr_id: " << chr_id << "（窗口范围: " << start_window << "~" << end_window << "）" << std::endl;

        // 启动新线程处理当前染色体
        threads.emplace_back(
            process_window_range,
            chr_id,
            start_window,
            end_window,
            window_size,
            read_length,
            full_chrom_seq,
            chr_bg_wb_ID
        );
    }

    // 3.3 等待最后一批线程完成
    for (auto& t : threads) {
        t.join();
    }

    std::cout << "所有染色体处理完成，总窗口数: " << g_window_info.size() << std::endl;
    return g_window_info;
}

std::vector<int> build_chr_bg_wb_ID() {
    std::vector<int> chr_bg_wb_ID(genome.size() + 1, 0);
    for (size_t i = 0; i < genome.size(); ++i) {
        const auto& chr_pair = genome[i];
        const std::string& current_id = chr_pair.first;
        const std::string& current_seq = chr_pair.second; 
        int seq_len = current_seq.size();
        int window_cnt = (seq_len - g_global_read_length) / g_global_read_length + 1;
        chr_bg_wb_ID[i+1] = chr_bg_wb_ID[i] + window_cnt + 5;
    }
    return chr_bg_wb_ID;
}

void BWT_CLASSIFY_MAIN::init_run(int argc, char *argv[]){
    double cpu_time = cputime();
	double real_time1 = realtime();
    share = (CLASSIFY_SHARE_DATA*)xcalloc(1, sizeof(CLASSIFY_SHARE_DATA));
	g_share = share;
    share->o = (OL_PAR *)xcalloc(1, sizeof(OL_PAR));
    if (share->o->get_option(argc, argv) != 0) {
        return;
    }
    g_input_is_paired = (share->o->paired_end != 0);
    char *input_bam_fn = share->o->read_bam;
    htsFile *input_file = hts_open(input_bam_fn, "rb");
    bam_hdr_t *header = sam_hdr_read(input_file);
    share->header = header;
    fprintf(stderr, "BEGIN LOADING INDEX\n");
    share->idx = (Aln_online::IDX_loader *)xcalloc(1, sizeof(Aln_online::IDX_loader));
    int read_length = share->o->read_length;
	genome = read_fasta(std::string(share->o->fasta_path));
    g_global_idx = share->idx;
    g_global_o = share->o;
    g_global_read_length = read_length;

    chr_bg_wb_ID = build_chr_bg_wb_ID();
    std::string work_dir = share->o->work_dir;
    int window_file_partition = chr_bg_wb_ID[genome.size()] / 10000 + 1;  // 分区数量
    share->blockFiles.clear();  // 确保文件句柄列表为空
    for (int i = 0; i <= window_file_partition; ++i) {
        const std::string file_name = work_dir + "/tmp/test." + std::to_string(i) + ".bin";
        FILE *fp = fopen(file_name.c_str(), "wb");  // 二进制写入模式
        if (!fp) {
            fprintf(stderr, "Failed to create partition file: %s\n", file_name.c_str());
            exit(EXIT_FAILURE);
        }
        share->blockFiles.emplace_back(fp);
    }
    
	g_partition_locks.clear();  // 确保为空
    for (size_t i = 0; i <=window_file_partition; ++i) {
        g_partition_locks.emplace_back(std::make_unique<std::mutex>());
    }
    fprintf(stderr, "Created %d partition files for window_id grouping\n", window_file_partition + 1);

    init_error_fastq_writer();

    // 4. 启动线程
    int thread_num = share->o->thread_n;
    std::vector<std::thread> consumer_threads;
    
	g_window_info.resize(chr_bg_wb_ID[genome.size()] + 1);
    g_window_info = generate_window_info(chr_bg_wb_ID, thread_num);  // 生成映射文件

    // 启动消费者线程（N个）
    for (int i = 0; i < thread_num; i++) {
        consumer_threads.emplace_back(consumer_thread);
    }

	std::thread producer(
        producer_thread,       // 你的生产者线程函数
        input_file,            // BAM文件句柄（生产者负责读取）
        header,                // BAM头信息
        share->o->work_dir     // 工作目录（传给BAMTask）
    );
    
    producer.join();

    // 9. 等待消费者线程结束（同原代码）
    for (auto& t : consumer_threads) {
        t.join();
    }

    close_error_fastq_writer();
    
	for (FILE *fp : share->blockFiles) {
        if (fp) {
            fclose(fp);
        }
    }
    share->blockFiles.clear();
    fprintf(stderr, "All partition files closed\n");
    
    
    std::vector<uint32_t> last_window_id(chr_bg_wb_ID[genome.size()]/10000 + 1, 0);

    int max_threads = thread_num;  // 复用配置的线程数
    std::vector<std::thread> threads;
    std::mutex mtx;  // 用于同步日志输出

    
    // 启动线程（控制并发数不超过max_threads）
    int active_threads = 0;
    for (int i = 0; i < window_file_partition; ++i) {
        // 若线程数达到上限，等待一个线程完成
        if (active_threads >= max_threads) {
            for (auto it = threads.begin(); it != threads.end(); ++it) {
                if (it->joinable()) {
                    it->join();
                    threads.erase(it);
                    active_threads--;
                    break;
                }
            }
        }

        // 启动线程处理第i个分区
        threads.emplace_back(
            process_partition,
            i,
            std::cref(work_dir),
            std::ref(last_window_id),
            std::ref(mtx)
        );
        active_threads++;
    }

    // 等待所有剩余线程完成
    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
    
	const int partition_count = chr_bg_wb_ID[genome.size()]/10000 + 1;
    std::vector<uint32_t> pre_window_id(partition_count, 0);
    std::vector<uint32_t> pre_window_id_record(partition_count, 0);
    uint32_t pre = 0;

    for (int i = 0; i < window_file_partition; i++) {
        if (last_window_id[i] == 0) continue;
        else {
            pre_window_id[i] = pre;
            pre = last_window_id[i];
        }
    }
    for (int i = 0; i <= window_file_partition; i++) {
        pre_window_id_record[i] = pre_window_id[i];
    }

    parallel_sort_and_write(
        pre_window_id,
        window_file_partition + 1,
        share->o->pos_dir,
        share->o->quality_score_dir,
        share->o->webp_dir,
        work_dir,
        share->o->thread_n
    );
	
    std::string parent_dir_path = get_parent_dir(share->o->pos_dir);
    if (share->o->webp_lossless) {
        std::cout << "[Quality] zstd lossless mode: quality_score directory will be packed into .xzip" << std::endl;
    } else {
        std::cout << "[Quality] WebP lossy mode: encode quality scores with WebP" << std::endl;
        webp_main_40(share->o->quality_score_dir, share->o->webp_dir, share->o->thread_n,
          chr_bg_wb_ID[genome.size()]/10000, g_global_read_length, share->o->webp_lossless, share->o->webp_quality,
          share->o->webp_method, share->o->webp_width, share->o->webp_height);
    }
    
    std::string pre_window_id_txt_path = parent_dir_path + "/pre_window_id.txt";
    
    std::ofstream pre_window_txt_file(pre_window_id_txt_path, std::ios::out | std::ios::trunc);
    if (pre_window_txt_file.is_open()) {
        int pre_window_size = window_file_partition + 1;
        for (int i = 0; i < pre_window_size; ++i) {
            pre_window_txt_file << pre_window_id_record[i] << std::endl;
        }
        pre_window_txt_file.close();
        std::cout << "pre_window_id 已保存到 " << pre_window_id_txt_path
                  << "，共 " << pre_window_size << " 行。" << std::endl;
    } else {
        std::cerr << "无法写入文本文件 " << pre_window_id_txt_path << std::endl;
    }
     
    std::vector<std::string> archive_paths = {
        "final_store_with_name",
        "diff_seq",
        "diff_base",
        "byte_flags",
        "pre_window_id.txt",
        "error.fastq"
    };
    if (!share->o->discard_qname) {
        archive_paths.push_back("qnames");
    }
    if (share->o->webp_lossless) {
        archive_paths.push_back("quality_score");
    } else {
        archive_paths.push_back("store_2_zero_lossless");
        archive_paths.push_back("quality_diff_base");
    }

    std::string xzip_path = parent_dir_path + "/compressed.xzip";
    if (!create_xzip_archive(parent_dir_path, xzip_path, archive_paths,
                             share->o->thread_n, share->o->zstd_level)) {
        std::cerr << "[XZIP] failed to create " << xzip_path << std::endl;
    }

    // 8. 释放所有资源（同原代码）
    bam_hdr_destroy(header);
    free(share->o);
    free(share->idx);
    free(share);
    
    
    
    fprintf(stderr, "Total CPU time: %.2f s\n", cputime() - cpu_time);
	fprintf(stderr, "Total REAL time: %.2f s\n", realtime() - real_time1);
}


std::atomic<uint32_t> g_pre_hap_id(0);

struct Task {
    int file_idx;
    const char* pos_dir;
    std::vector<uint32_t> pre_window_id_vector;
    CLASSIFY_SHARE_DATA* share;
    const std::vector<std::pair<std::string, std::string>>* genome;

    Task() = default;
    Task(int idx,
         const char* pd,
         std::vector<uint32_t> pre_wid,
         CLASSIFY_SHARE_DATA* s,
         const std::vector<std::pair<std::string, std::string>>* g)
        : file_idx(idx),
          pos_dir(pd),
          pre_window_id_vector(std::move(pre_wid)),
          share(s),
          genome(g) {}
};

// 存储单条read的完整信息
struct ReadPair {
    std::string qname;          // 序列名（来自qnames）
    std::string seq1;           // 第一条序列（来自reverse_seq）
    std::string seq2;           // 第二条序列（来自reverse_seq）
    std::string qual1;          // 第一条质量分数（来自quality_score）
    std::string qual2;          // 第二条质量分数（来自quality_score）
    uint16_t flag1;             // 第一条Flags（来自flags）
    uint16_t flag2;             // 第二条Flags（来自flags）
};

// 反向互补DNA序列（A<->T, C<->G，同时反转）
std::string reverse_complement(const std::string& seq) {
    std::string rev_comp;
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        switch (*it) {
            case 'A': rev_comp += 'T'; break;
            case 'T': rev_comp += 'A'; break;
            case 'C': rev_comp += 'G'; break;
            case 'G': rev_comp += 'C'; break;
            default: rev_comp += *it;  // 保留N等其他字符
        }
    }
    return rev_comp;
}

// 反向质量分数（仅反转，不互补）
std::string reverse_quality(const std::string& qual) {
    return std::string(qual.rbegin(), qual.rend());
}

std::vector<Compress_block> read_from_bin(const std::string& filename, uint32_t prev_window_id)
{
    Compress_final_block_with_huffman_table loaded;
    std::vector<Compress_bio_string_block> bio_string_blocks;
    std::vector<Compress_block> origin_blocks;
    std::vector<Compress_block> next_blocks;

    if (load_from_file(filename, loaded)) {
        std::vector<Compress_final_block> blocks = loaded.Final_blocks;
        std::unordered_map<char, std::string> huffmantable = loaded.huffmanCode;
        std::cout << "此处的blocks size为: " << blocks.size() << std::endl;

        std::unordered_map<std::string, char> reverse_map;
        for (const auto& pair : huffmantable) {
            reverse_map[pair.second] = pair.first;
        }

        for (auto& block : blocks) {
            Compress_bio_string_block bio_string_block;
            bio_string_block.bio_string = decode_huffman(block.huffman_code, reverse_map);
            memcpy(bio_string_block.qname, block.qname, block.lqname);
            bio_string_block.lqname    = block.lqname;
            bio_string_block.real_seq1 = block.real_seq1;
            bio_string_block.real_seq2 = block.real_seq2;
            bio_string_blocks.push_back(bio_string_block);
        }

        next_blocks = bio_string_to_num(bio_string_blocks, &origin_blocks, prev_window_id);
        std::cout << "文件 " << filename << " 处理完成" << std::endl;
    }
    return next_blocks;
}

std::vector<std::pair<std::string, std::string>> read_quality_scores(
    const std::string& base_dir,
    int id,
    const std::vector<Compress_block>& blocks,
    bool quality_lossless) {
    std::vector<std::pair<std::string, std::string>> quality_pairs;
    const std::string quality_dir = quality_lossless ? "quality_score" : "output_dir";
    std::string file_path = base_dir + "/" + quality_dir + "/quality_score." + std::to_string(id) + ".bin";
    std::ifstream in(file_path, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "无法打开质量分数文件: " << file_path << std::endl;
        return quality_pairs;
    }

    size_t block_count;
    in.read(reinterpret_cast<char*>(&block_count), sizeof(size_t));
    size_t expected_quality_count = 0;
    for (const auto& block : blocks) {
        expected_quality_count += block.is_paired ? 2 : 1;
    }
    if (block_count != expected_quality_count) {
        std::cerr << "Warning: quality count mismatch in " << file_path
                  << ", file=" << block_count
                  << ", expected=" << expected_quality_count << std::endl;
    }

    for (const auto& block : blocks) {
        // 读取第一条质量分数
        size_t len1;
        if (!in.read(reinterpret_cast<char*>(&len1), sizeof(size_t))) break;
        std::string qual1(len1, '\0');
        in.read(&qual1[0], len1);

        std::string qual2;
        if (block.is_paired) {
            // 读取第二条质量分数
            size_t len2;
            if (!in.read(reinterpret_cast<char*>(&len2), sizeof(size_t))) break;
            qual2.resize(len2);
            in.read(&qual2[0], len2);
        }
        for(int j=0;j<len1;j++)
        {
            if((int)qual1[j]<33)qual1[j]='!';
            if((int)qual1[j]>93)qual1[j]=']';
        }
        for(size_t j=0;j<qual2.size();j++)
        {
            if((int)qual2[j]<33)qual2[j]='!';
            if((int)qual2[j]>93)qual2[j]=']';
        }
        quality_pairs.emplace_back(qual1, qual2);
    }
    return quality_pairs;
}

// 从reverse_*.txt读取回复序列
std::vector<std::pair<std::string, std::string>> read_reverse_seqs(const std::string& base_dir, int id) {
    std::vector<std::pair<std::string, std::string>> seq_pairs;
    std::string file_path = base_dir + "/reverse_seq/reverse_" + std::to_string(id) + ".txt";
    std::ifstream in(file_path);
    if (!in.is_open()) {
        std::cerr << "无法打开回复序列文件: " << file_path << std::endl;
        return seq_pairs;
    }

    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        size_t index;
        std::string seq1, seq2;
        if (iss >> index >> seq1 >> seq2) {
            seq_pairs.emplace_back(seq1, seq2);
        }
    }
    return seq_pairs;
}

std::vector<std::pair<uint16_t, uint16_t>> read_flags(const std::string& parent_dir, int id)
{
    std::string file_path = parent_dir + "/byte_flags/flags_" + std::to_string(id) + ".bin";
    std::ifstream in(file_path, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "无法打开flags文件: " << file_path << std::endl;
        return {};
    }

    uint64_t pair_count = 0;
    in.read(reinterpret_cast<char*>(&pair_count), sizeof(pair_count));

    const size_t byte_count = (pair_count * 2 + 7) / 8;
    std::vector<uint8_t> packed(byte_count, 0);
    in.read(reinterpret_cast<char*>(packed.data()), byte_count);
    in.close();

    // 解包：只恢复 0x10 位，其余位全零（解压侧只用 & 0x10，完全够用）
    std::vector<std::pair<uint16_t, uint16_t>> flags;
    flags.reserve(pair_count);

    for (uint64_t i = 0; i < pair_count; ++i) {
        uint16_t f1 = (packed[(i * 2)     / 8] >> ((i * 2)     % 8)) & 1u ? 0x10 : 0x00;
        uint16_t f2 = (packed[(i * 2 + 1) / 8] >> ((i * 2 + 1) % 8)) & 1u ? 0x10 : 0x00;
        flags.emplace_back(f1, f2);
    }

    return flags;
}

// 从qnames_*.txt读取序列名
std::vector<std::string> read_qnames(const std::string& base_dir, int id) {
    std::vector<std::string> qnames;
    std::string file_path = base_dir + "/qnames/qnames_" + std::to_string(id) + ".txt";
    std::ifstream in(file_path);
    if (!in.is_open()) {
        std::cerr << "无法打开序列名文件: " << file_path << std::endl;
        return qnames;
    }

    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        size_t index;
        std::string qname;
        if (iss >> index >> qname) {
            qnames.push_back(qname);
        }
    }
    return qnames;
}

std::vector<std::string> generate_deterministic_qnames(int id, size_t count) {
    std::vector<std::string> qnames;
    qnames.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        qnames.push_back("XZIP_" + std::to_string(id) + "_" + std::to_string(i));
    }
    return qnames;
}

bool is_valid_fastq_4line(const std::string& fastq_str, const std::string& read_id, int pair_idx) {
    std::istringstream ss(fastq_str);
    std::string line;
    int line_count = 0;
    
    // 统计有效行数（排除空行，避免换行符问题）
    while (std::getline(ss, line)) {
        // 去除行尾可能的\r（Windows换行符）
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        // 空行直接判定为格式错误
        if (line.empty()) {
            std::cerr << "[ERROR] 第" << pair_idx << "对read的" << read_id 
                      << "包含空行，不符合FASTQ格式" << std::endl;
            return false;
        }
        line_count++;
        
        // 按FASTQ规范校验每行格式（给每个case加{}作用域）
        switch (line_count) {
            case 1: {  // 加{}
                if (line.empty() || line[0] != '@') {
                    std::cerr << "[ERROR] 第" << pair_idx << "对read的" << read_id 
                              << "首行未以@开头: " << line << std::endl;
                    return false;
                }
                break;
            }
            case 2: {  // 加{}
                // 序列行不能包含非法字符（仅允许ATCGN/atcgn）
                for (char c : line) {
                    if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == 'N' || c == 'n')) {
                        std::cerr << "[ERROR] 第" << pair_idx << "对read的" << read_id 
                                  << "序列行包含非法字符: " << c << std::endl;
                        return false;
                    }
                }
                break;
            }
            case 3: {  // 加{}
                if (line != "+") {
                    std::cerr << "[ERROR] 第" << pair_idx << "对read的" << read_id 
                              << "第三行必须是+: " << line << std::endl;
                    return false;
                }
                break;
            }
            case 4: {  // 加{}（关键：这里定义了局部变量，必须加{}）
                // 质量行长度需与序列行一致
                std::istringstream ss_check(fastq_str);
                std::string seq_line;
                std::getline(ss_check, seq_line); // 跳过首行
                std::getline(ss_check, seq_line); // 读取序列行
                if (line.length() != seq_line.length()) {
                    std::cerr << "[ERROR] 第" << pair_idx << "对read的" << read_id 
                              << "质量行长度(" << line.length() << ")与序列行长度(" 
                              << seq_line.length() << ")不匹配" << std::endl;
                    return false;
                }
                break;
            }
            default: {  // 加{}
                std::cerr << "[ERROR] 第" << pair_idx << "对read的" << read_id 
                          << "行数超过4行（实际" << line_count << "行）" << std::endl;
                return false;
            }
        }
    }

    // 校验行数是否为4行
    if (line_count != 4) {
        std::cerr << "[ERROR] 第" << pair_idx << "对read的" << read_id 
                  << "行数不符合要求（期望4行，实际" << line_count << "行）" << std::endl;
        return false;
    }

    return true;
}

void append_error_fastq_interleaved(
    const std::string& error_fastq_path,
    std::ofstream& merged_r1,
    std::ofstream& merged_r2
) {
    std::ifstream in(error_fastq_path);
    if (!in.is_open()) {
        std::cerr << "[WARN] no error.fastq found: " << error_fastq_path << std::endl;
        return;
    }

    std::string h1, s1, p1, q1;
    std::string h2, s2, p2, q2;
    uint64_t error_pairs = 0;

    while (std::getline(in, h1)) {
        if (!std::getline(in, s1) ||
            !std::getline(in, p1) ||
            !std::getline(in, q1)) {
            std::cerr << "[WARN] truncated R1 record in error.fastq" << std::endl;
            break;
        }

        if (!std::getline(in, h2) ||
            !std::getline(in, s2) ||
            !std::getline(in, p2) ||
            !std::getline(in, q2)) {
            std::cerr << "[WARN] truncated R2 record in error.fastq" << std::endl;
            break;
        }

        merged_r1 << h1 << "\n" << s1 << "\n" << p1 << "\n" << q1 << "\n";
        merged_r2 << h2 << "\n" << s2 << "\n" << p2 << "\n" << q2 << "\n";

        ++error_pairs;
    }

    std::cout << "[MERGE] appended error.fastq pairs = "
              << error_pairs << std::endl;
}

void append_error_fastq_single(
    const std::string& error_fastq_path,
    std::ofstream& merged_r1
) {
    std::ifstream in(error_fastq_path);
    if (!in.is_open()) {
        std::cerr << "[WARN] no error.fastq found: " << error_fastq_path << std::endl;
        return;
    }

    std::string h, s, p, q;
    uint64_t error_reads = 0;
    while (std::getline(in, h)) {
        if (!std::getline(in, s) ||
            !std::getline(in, p) ||
            !std::getline(in, q)) {
            std::cerr << "[WARN] truncated single-end record in error.fastq" << std::endl;
            break;
        }
        merged_r1 << h << "\n" << s << "\n" << p << "\n" << q << "\n";
        ++error_reads;
    }

    std::cout << "[MERGE] appended error.fastq reads = "
              << error_reads << std::endl;
}

void go_to_decompress_pos(std::vector<uint32_t> pre_window_id_vector,
    CLASSIFY_SHARE_DATA* share)
{
    std::atomic<uint64_t> block_count_sum(0);
    std::atomic<uint64_t> correct(0);
    double start_time = cputime();

    int read_length = share->o->read_length;

    constexpr bool ONLY_RUN_TARGET_TASK = false;
    constexpr int  TARGET_TASK_IDX = 0;
    constexpr bool MERGE_ONLY_TARGET_TASK = true;
    constexpr bool DECOMPRESS_DEBUG_SINGLE_THREAD = false;
    constexpr bool SIMPLE_BLOCK_PROGRESS = true;
    constexpr int  BLOCK_LOG_INTERVAL = 10000;

    std::queue<Task> task_queue;
    std::mutex queue_mutex;
    std::condition_variable cv;
    std::atomic<bool> stop_flag(false);
    std::atomic<int> completed_tasks(0);

    int file_num = chr_bg_wb_ID[genome.size()] / 10000;

    std::vector<int> task_ids_to_run;
    if (ONLY_RUN_TARGET_TASK) {
        if (TARGET_TASK_IDX >= 0 && TARGET_TASK_IDX <= file_num) {
            task_ids_to_run.push_back(TARGET_TASK_IDX);
        } else {
            std::cerr << "[FATAL] TARGET_TASK_IDX=" << TARGET_TASK_IDX
                      << " 超出范围 [0, " << file_num << "]，函数直接返回。" << std::endl;
            return;
        }
    } else {
        for (int i = 0; i <= file_num; ++i) {
            task_ids_to_run.push_back(i);
        }
    }

    int total_tasks = static_cast<int>(task_ids_to_run.size());
    if (total_tasks == 0) {
        std::cerr << "[WARN] 没有可运行的 task，直接返回。" << std::endl;
        return;
    }

    std::cout << "[MAIN] file_num=" << file_num
              << ", total_tasks=" << total_tasks
              << ", read_length=" << read_length
              << ", ONLY_RUN_TARGET_TASK=" << ONLY_RUN_TARGET_TASK
              << std::endl;

    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        for (int idx : task_ids_to_run) {
            task_queue.push(Task(
                idx,
                share->o->pos_dir,
                pre_window_id_vector,
                share,
                &genome
            ));
        }
    }

    std::cout << "[MAIN] queued tasks = " << task_ids_to_run.size() << std::endl;

    auto process_task = [&]() {
        // 通过 window_id 切参考片段
        auto get_ref_by_window_id = [&](uint32_t window_id, std::string& ref_out) -> bool {
            if (window_id >= g_window_info.size()) return false;

            const Window_t& w = g_window_info[window_id];
            int chr_ID = w.chr_ID;
            uint32_t st_pos = w.st_pos;

            if (chr_ID < 0 || chr_ID >= (int)genome.size()) return false;
            if (st_pos == 0) return false;

            const std::string& full_chrom_seq = genome[chr_ID].second;
            uint32_t start_idx = st_pos - 1;
            if (start_idx >= full_chrom_seq.size()) return false;

            int actual_length = std::min(2 * read_length, (int)(full_chrom_seq.size() - start_idx));
            if (actual_length <= 0) return false;

            ref_out.assign(full_chrom_seq, start_idx, actual_length);
            return true;
        };

        // 重建单条 read
        auto reconstruct_read = [&](
            const std::string& diff_seq,
            const std::string& diff_base_all,
            int& diff_base_idx,
            uint32_t window_id,
            uint32_t hap_offset,
            std::string& out_read,
            const std::string& task_label,
            int block_idx,
            bool& fatal_error
        ) {
            std::string ref;
            if (!get_ref_by_window_id(window_id, ref)) {
                fprintf(stderr, "[ERROR] %s get_ref failed block_idx=%d\n",
                        task_label.c_str(), block_idx);
                fatal_error = true;
                return;
            }

            if ((int)diff_seq.size() < read_length) {
                fprintf(stderr, "[ERROR] %s diff_seq.size=%zu < read_length=%d block_idx=%d\n",
                        task_label.c_str(), diff_seq.size(), read_length, block_idx);
                fatal_error = true;
                return;
            }

            out_read.reserve(read_length);
            for (int k = 0; k < read_length; ++k) {
                if (diff_seq[k] == '1') {
                    if (diff_base_idx >= (int)diff_base_all.size()) {
                        fprintf(stderr, "[ERROR] %s diff_base 越界 block_idx=%d k=%d\n",
                                task_label.c_str(), block_idx, k);
                        fatal_error = true;
                        return;
                    }
                    out_read += diff_base_all[diff_base_idx++];
                } else {
                    int ref_k = (int)hap_offset + k;
                    if (ref_k >= (int)ref.size()) {
                        fprintf(stderr, "[ERROR] %s ref 越界 block_idx=%d ref_k=%d ref.size=%zu\n",
                                task_label.c_str(), block_idx, ref_k, ref.size());
                        fatal_error = true;
                        return;
                    }
                    out_read += ref[ref_k];
                }
            }
        };

        auto err = [&](const std::string& msg, int task_file_idx) {
            std::cerr << "[ERROR][thread=" << std::this_thread::get_id()
                      << "][task=" << task_file_idx << "] "
                      << msg << std::endl;
        };

        while (true) {
            Task task;
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                cv.wait(lock, [&]() {
                    return stop_flag || !task_queue.empty();
                });

                if (stop_flag && task_queue.empty()) return;
                if (task_queue.empty()) continue;

                task = std::move(task_queue.front());
                task_queue.pop();
            }

            std::cout << "[TASK-START] thread=" << std::this_thread::get_id()
                      << " task_idx=" << task.file_idx
                      << std::endl;

            std::string filename = std::string(task.pos_dir) + "/store." +
                                   std::to_string(task.file_idx) + ".bin";

            std::vector<Compress_block> blocks = read_from_bin(
                filename,
                pre_window_id_vector[task.file_idx]);

            std::string parent_dir = get_parent_dir(task.pos_dir);

            std::cout << "[LOAD] task_idx=" << task.file_idx
                      << " blocks.size()=" << blocks.size()
                      << std::endl;

            std::string diff_base_decoded_str = diff_seq_read(task.pos_dir, task.file_idx);
            if (diff_base_decoded_str.empty()) {
                err("diff_base 解码失败或为空。", task.file_idx);
                int done = ++completed_tasks;
                std::cout << "[TASK-END] task_idx=" << task.file_idx
                          << " status=diff_base_empty completed_tasks=" << done << "/" << total_tasks
                          << std::endl;
                cv.notify_all();
                continue;
            }

            std::cout << "[LOAD] task_idx=" << task.file_idx
                      << " diff_base.size()=" << diff_base_decoded_str.size()
                      << std::endl;

            std::vector<std::string> diff_seqs = read_diff_seq_blocks(task.pos_dir, task.file_idx);
            if (diff_seqs.empty()) {
                err("diff_seq 读取失败或为空。", task.file_idx);
                int done = ++completed_tasks;
                std::cout << "[TASK-END] task_idx=" << task.file_idx
                          << " status=diff_seq_empty completed_tasks=" << done << "/" << total_tasks
                          << std::endl;
                cv.notify_all();
                continue;
            }

            std::cout << "[LOAD] task_idx=" << task.file_idx
                      << " diff_seqs.size()=" << diff_seqs.size()
                      << std::endl;

            if (blocks.empty()) {
                int done = ++completed_tasks;
                std::cout << "[TASK-END] task_idx=" << task.file_idx
                          << " status=blocks_empty completed_tasks=" << done << "/" << total_tasks
                          << std::endl;
                cv.notify_all();
                continue;
            }

            std::vector<std::pair<std::string, std::string>> reverse_pairs;
            reverse_pairs.reserve(blocks.size());

            int diff_base_idx = 0;
            int process_idx   = 0;
            int block_idx     = 0;
            bool fatal_error  = false;

            for (Compress_block& block : blocks) {
                ++block_idx;
                std::string reverse_read, reverse_read2;
                const bool block_is_paired = block.is_paired != 0;
                

                // ---------- R1 ----------
                if (process_idx >= (int)diff_seqs.size()) {
                    err("R1: process_idx=" + std::to_string(process_idx) +
                        " 超出 diff_seqs 范围 " + std::to_string(diff_seqs.size()), task.file_idx);
                    fatal_error = true;
                    break;
                }
                if(strcmp(block.qname,"A00900:594:HJVHNDSX7:1:2417:24840:13667")==0)
                {
                    std::cout<<"debug\n";
                }
                reconstruct_read(diff_seqs[process_idx], diff_base_decoded_str, diff_base_idx,
                                 block.window_id, block.hap_offset,
                                 reverse_read, "R1", block_idx, fatal_error);
                if (fatal_error) break;
                process_idx++;

                if (block_is_paired) {
                    // ---------- R2 ----------
                    if (process_idx >= (int)diff_seqs.size()) {
                        err("R2: process_idx=" + std::to_string(process_idx) +
                            " 超出 diff_seqs 范围 " + std::to_string(diff_seqs.size()), task.file_idx);
                        fatal_error = true;
                        break;
                    }
                    if(strcmp(block.qname,"A00900:594:HJVHNDSX7:1:2417:24840:13667")==0)
                    {
                        std::cout<<"debug\n";
                    }
                    reconstruct_read(diff_seqs[process_idx], diff_base_decoded_str, diff_base_idx,
                                     block.window_id2, block.hap_offset2,
                                     reverse_read2, "R2", block_idx, fatal_error);
                    if (fatal_error) break;
                    process_idx++;
                }
                /*
                ////debug正确性开始
                if (reverse_read != block.real_seq1) {
                    
                    fprintf(stderr, "[MISMATCH-R1] block_idx=%d\n", block_idx);
                    fprintf(stderr, "block.qname=%s\n", block.qname);
                    fprintf(stderr, "  real_seq1: %s\n", block.real_seq1.c_str());
                    fprintf(stderr, "  recovered: %s\n", reverse_read.c_str());
                    fprintf(stderr, "  window_id=%u hap_offset=%u\n", block.window_id, block.hap_offset);
                    fprintf(stderr, "  diff_seq1: %s\n", diff_seqs[process_idx-2].c_str());
                    // 逐位找第一个不同的位置
                    for (int ki = 0; ki < (int)std::min(block.real_seq1.size(), reverse_read.size()); ++ki) {
                        if (block.real_seq1[ki] != reverse_read[ki]) {
                            fprintf(stderr, "  first diff at k=%d real=%c recovered=%c\n",
                                    ki, block.real_seq1[ki], reverse_read[ki]);
                            break;
                        }
                    }
                }
                if (reverse_read2 != block.real_seq2) {
                    fprintf(stderr, "[MISMATCH-R2] block_idx=%d\n", block_idx);
                    fprintf(stderr, "block.qname=%s\n", block.qname);
                    fprintf(stderr, "  real_seq2: %s\n", block.real_seq2.c_str());
                    fprintf(stderr, "  recovered: %s\n", reverse_read2.c_str());
                    fprintf(stderr, "  window_id2=%u hap_offset2=%u\n", block.window_id2, block.hap_offset2);
                    fprintf(stderr, "  diff_seq2: %s\n", diff_seqs[process_idx-1].c_str());
                    for (int ki = 0; ki < (int)std::min(block.real_seq2.size(), reverse_read2.size()); ++ki) {
                        if (block.real_seq2[ki] != reverse_read2[ki]) {
                            fprintf(stderr, "  first diff at k=%d real=%c recovered=%c\n",
                                    ki, block.real_seq2[ki], reverse_read2[ki]);
                            break;
                        }
                    }
                }
                ////debug正确性结束
                */
                reverse_pairs.emplace_back(std::move(reverse_read), std::move(reverse_read2));

                if (SIMPLE_BLOCK_PROGRESS &&
                    (block_idx % BLOCK_LOG_INTERVAL == 0 || block_idx == (int)blocks.size())) {
                    std::cout << "[BLOCK] thread=" << std::this_thread::get_id()
                              << " task_idx=" << task.file_idx
                              << " block_idx=" << block_idx << "/" << blocks.size()
                              << " process_idx=" << process_idx
                              << " diff_base_idx=" << diff_base_idx
                              << std::endl;
                }
            }

            if (fatal_error) {
                err("因致命错误跳过，已处理 " + std::to_string(reverse_pairs.size()) + " 对。", task.file_idx);
                int done = ++completed_tasks;
                std::cout << "[TASK-END] task_idx=" << task.file_idx
                          << " status=fatal_error completed_tasks=" << done << "/" << total_tasks
                          << std::endl;
                cv.notify_all();
                continue;
            }

            auto quality_pairs = read_quality_scores(parent_dir, task.file_idx, blocks, task.share->o->webp_lossless != 0);
            auto qnames = task.share->o->discard_qname
                        ? generate_deterministic_qnames(task.file_idx, reverse_pairs.size())
                        : read_qnames(parent_dir, task.file_idx);
            auto flag_pairs    = read_flags(parent_dir, task.file_idx);
            bool task_has_paired_reads = false;
            for (const auto& block : blocks) {
                if (block.is_paired) {
                    task_has_paired_reads = true;
                    break;
                }
            }

            size_t total_pairs = reverse_pairs.size();
            bool meta_ok = (qnames.size()        == total_pairs &&
                            quality_pairs.size() == total_pairs &&
                            flag_pairs.size()    == total_pairs);

            if (!meta_ok) {
                err("数据量不匹配，reverse=" + std::to_string(total_pairs) +
                    ", quality=" + std::to_string(quality_pairs.size()) +
                    ", qnames=" + std::to_string(qnames.size()) +
                    ", flags=" + std::to_string(flag_pairs.size()), task.file_idx);
                int done = ++completed_tasks;
                std::cout << "[TASK-END] task_idx=" << task.file_idx
                          << " status=meta_mismatch completed_tasks=" << done << "/" << total_tasks
                          << std::endl;
                cv.notify_all();
                continue;
            }

            std::string final_dir = parent_dir + "/final_fastq";
            if (!create_directory(final_dir)) {
                err("无法创建 final_fastq 文件夹。", task.file_idx);
                int done = ++completed_tasks;
                std::cout << "[TASK-END] task_idx=" << task.file_idx
                          << " status=mkdir_failed completed_tasks=" << done << "/" << total_tasks
                          << std::endl;
                cv.notify_all();
                continue;
            }

            std::string fastq_R1_path = final_dir + "/output_" + std::to_string(task.file_idx) + "_R1.fastq";
            std::string fastq_R2_path = final_dir + "/output_" + std::to_string(task.file_idx) + "_R2.fastq";
            std::ofstream out_r1(fastq_R1_path);
            std::ofstream out_r2;
            if (task_has_paired_reads) {
                out_r2.open(fastq_R2_path);
            }

            if (!out_r1.is_open() || (task_has_paired_reads && !out_r2.is_open())) {
                err("无法创建 FASTQ 文件: " + fastq_R1_path + " 或 " + fastq_R2_path, task.file_idx);
                int done = ++completed_tasks;
                std::cout << "[TASK-END] task_idx=" << task.file_idx
                          << " status=fastq_open_failed completed_tasks=" << done << "/" << total_tasks
                          << std::endl;
                cv.notify_all();
                continue;
            }

            std::cout << "[WRITE] task_idx=" << task.file_idx
                      << " total_pairs=" << total_pairs
                      << " begin write fastq"
                      << std::endl;

            int error_count = 0;
            const int MAX_ERRORS = 100;

            for (size_t i = 0; i < total_pairs; ++i) {
                ReadPair rp;
                rp.qname = qnames[i];
                rp.seq1  = reverse_pairs[i].first;
                rp.seq2  = reverse_pairs[i].second;
                rp.qual1 = quality_pairs[i].first;
                rp.qual2 = quality_pairs[i].second;
                rp.flag1 = flag_pairs[i].first;
                rp.flag2 = flag_pairs[i].second;
                const bool read_is_paired = blocks[i].is_paired != 0;

                if (rp.flag1 & 0x10) {
                    rp.seq1  = reverse_complement(rp.seq1);
                    rp.qual1 = reverse_quality(rp.qual1);
                }
                if (read_is_paired && (rp.flag2 & 0x10)) {
                    rp.seq2  = reverse_complement(rp.seq2);
                    rp.qual2 = reverse_quality(rp.qual2);
                }

                std::string read1 = "@" + rp.qname + (read_is_paired ? "/1\n" : "\n") + rp.seq1 + "\n+\n" + rp.qual1 + "\n";
                std::string read2 = "@" + rp.qname + "/2\n" + rp.seq2 + "\n+\n" + rp.qual2 + "\n";

                if (!is_valid_fastq_4line(read1, "R1", i) ||
                    (read_is_paired && !is_valid_fastq_4line(read2, "R2", i))) {
                    if (++error_count >= MAX_ERRORS) {
                        err("错误数超过阈值，停止该 file_idx。", task.file_idx);
                        break;
                    }
                    continue;
                }

                out_r1 << read1;
                if (read_is_paired) out_r2 << read2;
            }

            out_r1.close();
            if (task_has_paired_reads) out_r2.close();

            std::cout << "[WRITE] task_idx=" << task.file_idx
                      << " finish write fastq"
                      << std::endl;

            int done = ++completed_tasks;
            std::cout << "[TASK-END] task_idx=" << task.file_idx
                      << " status=success completed_tasks=" << done << "/" << total_tasks
                      << std::endl;
            cv.notify_all();
        }
    };

    unsigned int requested_threads = DECOMPRESS_DEBUG_SINGLE_THREAD
        ? 1u
        : (share->o->thread_n > 0 ? (unsigned int)share->o->thread_n : 1u);

    unsigned int thread_count = std::max(1u, std::min(requested_threads, (unsigned int)total_tasks));

    std::cout << "[MAIN] requested_threads=" << requested_threads
              << ", thread_count=" << thread_count
              << ", DECOMPRESS_DEBUG_SINGLE_THREAD=" << DECOMPRESS_DEBUG_SINGLE_THREAD
              << std::endl;

    std::vector<std::thread> thread_pool;
    thread_pool.reserve(thread_count);
    for (unsigned int i = 0; i < thread_count; ++i) {
        thread_pool.emplace_back(process_task);
    }
    cv.notify_all();

    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        cv.wait(lock, [&]() { return completed_tasks == total_tasks; });
        stop_flag = true;
    }
    cv.notify_all();

    for (auto& t : thread_pool) t.join();

    std::string parent_dir_final = get_parent_dir(share->o->pos_dir);
    std::string final_dir        = parent_dir_final + "/final_fastq";
    std::string merged_R1_path   = parent_dir_final + "/merged_all_R1.fastq";
    std::string merged_R2_path   = parent_dir_final + "/merged_all_R2.fastq";

    std::ofstream merged_r1(merged_R1_path, std::ios::trunc);
    std::ofstream merged_r2;
    if (g_input_is_paired) {
        merged_r2.open(merged_R2_path, std::ios::trunc);
    }

    if (!merged_r1.is_open() || (g_input_is_paired && !merged_r2.is_open())) {
        std::cerr << "[FATAL] 无法创建 merge 输出文件。" << std::endl;
        return;
    }

    auto append = [](std::ofstream& dst, const std::string& path) {
        std::ifstream src(path);
        if (src.is_open()) dst << src.rdbuf();
    };

    std::cout << "[MERGE] begin merge fastq" << std::endl;

    if (ONLY_RUN_TARGET_TASK && MERGE_ONLY_TARGET_TASK) {
        append(merged_r1, final_dir + "/output_" + std::to_string(TARGET_TASK_IDX) + "_R1.fastq");
        if (g_input_is_paired) {
            append(merged_r2, final_dir + "/output_" + std::to_string(TARGET_TASK_IDX) + "_R2.fastq");
        }
    } else {
        for (int id = 0; id <= file_num; ++id) {
            append(merged_r1, final_dir + "/output_" + std::to_string(id) + "_R1.fastq");
            if (g_input_is_paired) {
                append(merged_r2, final_dir + "/output_" + std::to_string(id) + "_R2.fastq");
            }
        }
    }

    std::string error_fastq_path = parent_dir_final + "/error.fastq";
    if (g_input_is_paired) {
        append_error_fastq_interleaved(error_fastq_path, merged_r1, merged_r2);
    } else {
        append_error_fastq_single(error_fastq_path, merged_r1);
    }

    merged_r1.close();
    if (g_input_is_paired) merged_r2.close();

    std::cout << "[MERGE] finish merge fastq" << std::endl;
    std::cout << "[MAIN] go_to_decompress_pos finished, elapsed="
              << (cputime() - start_time) << " sec" << std::endl;
}

void DECOMPRESS_MAIN::run(int argc, char* argv[])
{
    std::cout << "work decompress" << std::endl;
    double cpu_time = cputime();
    CLASSIFY_SHARE_DATA* share = NULL;
    share = (CLASSIFY_SHARE_DATA*)xcalloc(1, sizeof(CLASSIFY_SHARE_DATA));
    g_share = share;
    share->o = (OL_PAR*)xcalloc(1, sizeof(OL_PAR));
    if (share->o->get_decompress_option(argc, argv) != 0)
        return;
    g_input_is_paired = (share->o->paired_end != 0);
    char output_dir[512];
    if (create_output_dir_c(share->o->webp_dir, output_dir, sizeof(output_dir)) == 0) {
        printf("后续操作可用路径：%s\n", output_dir);
    }
    g_global_idx = share->idx;
    g_global_o = share->o;
    g_global_read_length = share->o->read_length;
    std::vector<uint32_t> pre_window_id_vector;  // 改名
    std::string parent_dir = get_parent_dir(share->o->pos_dir);
    std::string xzip_path = parent_dir + "/compressed.xzip";
    if ((!directory_exists(share->o->pos_dir) ||
         !file_exists(parent_dir + "/pre_window_id.txt")) &&
        file_exists(xzip_path)) {
        std::cout << "[XZIP] restoring compressed archive: " << xzip_path << std::endl;
        if (!restore_xzip_archive(parent_dir, xzip_path)) {
            std::cerr << "[FATAL] cannot restore " << xzip_path << std::endl;
            free(share->idx); free(share->o); free(share);
            return;
        }
    }

    std::string pre_window_id_txt_path = parent_dir + "/pre_window_id.txt";
    std::ifstream pre_window_txt_file(pre_window_id_txt_path);
    if (!pre_window_txt_file.is_open()) {
        std::cerr << "错误：无法读取 " << pre_window_id_txt_path << std::endl;
        free(share->idx); free(share->o); free(share);
        return;
    }
    std::string line;
    while (std::getline(pre_window_txt_file, line)) {
        try {
            pre_window_id_vector.push_back(static_cast<uint32_t>(std::stoul(line)));
        } catch (...) {}
    }
    pre_window_txt_file.close();
    if (pre_window_id_vector.empty()) {
        std::cerr << "错误：pre_window_id 为空" << std::endl;
        free(share->idx); free(share->o); free(share);
        return;
    }
    std::cout << "成功读取 " << pre_window_id_vector.size() << " 个 pre_window_id。" << std::endl;
    genome = read_fasta(std::string(share->o->fasta_path));
    chr_bg_wb_ID = build_chr_bg_wb_ID();
    g_window_info.resize(chr_bg_wb_ID[genome.size()] + 1);
    g_window_info = generate_window_info(chr_bg_wb_ID, share->o->thread_n);
    if (share->o->webp_lossless) {
        std::cout << "[Quality] zstd lossless mode: use quality_score directory from .xzip/work_dir" << std::endl;
        if (!directory_exists(parent_dir + "/quality_score") &&
            file_exists(xzip_path) &&
            !restore_xzip_archive(parent_dir, xzip_path)) {
            std::cerr << "[FATAL] 无法恢复 compressed.xzip" << std::endl;
            free(share->idx); free(share->o); free(share);
            return;
        }
    } else {
        std::cout << "[Quality] WebP lossy mode: reconstruct quality scores" << std::endl;
        webp_reconstructor_40_main(
                share->o->webp_dir,
                output_dir,
                share->o->thread_n,
                share->o->read_length,
                chr_bg_wb_ID[genome.size()] / 10000
        );
    }
    go_to_decompress_pos(pre_window_id_vector, share);
    std::cout << "run decompress successfully" << std::endl;
    free(share->idx);
    free(share->o);
    free(share);
}
