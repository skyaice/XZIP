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
 #include <sys/resource.h>
 #include <iostream>
 #include <sstream>
 #include <queue>
 #include <array>
 #include <unordered_map>
 #include <algorithm>
 #include <inttypes.h>
 #include <filesystem>
 #include <system_error>
#include <chrono>
#include <ctime>
#include <thread>
 #include "BWT_aln.hpp"
 #include "CPPLIB/tools.hpp"
 extern "C"{
 #include "clib/kthread.h"
 #include "clib/utils.h"
 #include "clib/desc.h"
 }
 #include "occ_RST_DEF.hpp"
 #include "./BWT_idx/var_map.hpp"
 #include "./webp_tool/compress_40_quality_score.h"
 #include "./webp_tool/webp_reconstructor_40.h"


 #include "./htslib/htslib/sam.h"

 using namespace BWT_aln;
 using FastaData = std::vector<std::pair<std::string, std::string>>;
 namespace fs = std::filesystem;

using PerfClock = std::chrono::steady_clock;
extern OL_PAR* g_global_o;

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
    double hap_ms = 0.0;          // read-block transformation
    double block_ms = 0.0;        // 填充 Compress_block + push local_blocks
    double flush_ms = 0.0;        // flush_local_blocks
    double cleanup_ms = 0.0;      // BAM 对象归还对象池等
    double total_batch_ms = 0.0;  // 整个 batch 总耗时
};

static std::mutex g_perf_print_mutex;


 int g_global_read_length = 0;
 std::vector<int>chr_bg_wb_ID;
 std::vector<Window_t> g_window_info;
std::array<std::atomic<bool>, 256> g_quality_score_seen{};

 FastaData genome;

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

BamReadMetadata capture_bam_metadata(const bam1_t* read) {
    BamReadMetadata metadata;
    if (!read) return metadata;
    metadata.core = read->core;
    metadata.quality_missing = read->core.l_qseq > 0 ? 1 : 0;
    const uint8_t* quality = bam_get_qual(read);
    for (int i = 0; i < read->core.l_qseq; ++i) {
        if (quality[i] != 0xff) {
            metadata.quality_missing = 0;
            break;
        }
    }
    if (read->core.n_cigar > 0) {
        const uint32_t* cigar = bam_get_cigar(read);
        metadata.cigar.assign(cigar, cigar + read->core.n_cigar);
    }
    const int aux_length = bam_get_l_aux(read);
    const uint8_t* aux = bam_get_aux(read);
    if (aux_length > 0) metadata.aux.assign(aux, aux + aux_length);
    return metadata;
}

void append_bam_metadata(std::vector<char>& output, const BamReadMetadata& metadata) {
    const char* core_bytes = reinterpret_cast<const char*>(&metadata.core);
    output.insert(output.end(), core_bytes, core_bytes + sizeof(metadata.core));
    output.push_back(static_cast<char>(metadata.quality_missing));

    const uint32_t cigar_count = static_cast<uint32_t>(metadata.cigar.size());
    const uint32_t aux_size = static_cast<uint32_t>(metadata.aux.size());
    output.insert(output.end(), reinterpret_cast<const char*>(&cigar_count),
                  reinterpret_cast<const char*>(&cigar_count) + sizeof(cigar_count));
    if (cigar_count > 0) {
        const char* bytes = reinterpret_cast<const char*>(metadata.cigar.data());
        output.insert(output.end(), bytes, bytes + cigar_count * sizeof(uint32_t));
    }
    output.insert(output.end(), reinterpret_cast<const char*>(&aux_size),
                  reinterpret_cast<const char*>(&aux_size) + sizeof(aux_size));
    if (aux_size > 0) {
        const char* bytes = reinterpret_cast<const char*>(metadata.aux.data());
        output.insert(output.end(), bytes, bytes + aux_size);
    }
}

void append_bam_record_metadata(std::vector<char>& output, const bam1_t* read) {
    if (!read) return;
    const char* core_bytes = reinterpret_cast<const char*>(&read->core);
    output.insert(output.end(), core_bytes, core_bytes + sizeof(read->core));
    uint8_t quality_missing = read->core.l_qseq > 0 ? 1 : 0;
    const uint8_t* quality = bam_get_qual(read);
    for (int i = 0; i < read->core.l_qseq; ++i) {
        if (quality[i] != 0xff) {
            quality_missing = 0;
            break;
        }
    }
    output.push_back(static_cast<char>(quality_missing));

    const uint32_t cigar_count = read->core.n_cigar;
    output.insert(output.end(), reinterpret_cast<const char*>(&cigar_count),
                  reinterpret_cast<const char*>(&cigar_count) + sizeof(cigar_count));
    if (cigar_count > 0) {
        const char* bytes = reinterpret_cast<const char*>(bam_get_cigar(read));
        output.insert(output.end(), bytes, bytes + cigar_count * sizeof(uint32_t));
    }

    const uint32_t aux_size = bam_get_l_aux(read);
    output.insert(output.end(), reinterpret_cast<const char*>(&aux_size),
                  reinterpret_cast<const char*>(&aux_size) + sizeof(aux_size));
    if (aux_size > 0) {
        const char* bytes = reinterpret_cast<const char*>(bam_get_aux(read));
        output.insert(output.end(), bytes, bytes + aux_size);
    }
}

bool write_bam_metadata_file(const std::string& path,
                             const std::vector<Compress_block>& blocks) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) return false;

    static const char magic[8] = {'X', 'B', 'A', 'M', 'M', 'E', 'T', '1'};
    const uint32_t version = 2;
    const uint64_t block_count = blocks.size();
    out.write(magic, sizeof(magic));
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    out.write(reinterpret_cast<const char*>(&block_count), sizeof(block_count));

    std::vector<char> buffer;
    buffer.reserve(4 * 1024 * 1024);
    for (const Compress_block& block : blocks) {
        buffer.push_back(static_cast<char>(block.is_paired));
        append_bam_metadata(buffer, block.bam_meta1);
        if (block.is_paired) append_bam_metadata(buffer, block.bam_meta2);
        if (buffer.size() >= 4 * 1024 * 1024) {
            out.write(buffer.data(), buffer.size());
            buffer.clear();
        }
    }
    if (!buffer.empty()) out.write(buffer.data(), buffer.size());
    return out.good();
}

bool file_exists(const std::string& file_path) {
    std::error_code ec;
    return fs::is_regular_file(file_path, ec);
}

bool is_zstd_frame(const std::string& file_path) {
    std::ifstream in(file_path, std::ios::binary);
    unsigned char magic[4] = {};
    if (!in.read(reinterpret_cast<char*>(magic), sizeof(magic))) return false;
    return magic[0] == 0x28 && magic[1] == 0xb5 &&
           magic[2] == 0x2f && magic[3] == 0xfd;
}

bool compress_bam_info_metadata(const std::string& bam_info_dir, int thread_n) {
    const int zstd_threads = std::max(1, thread_n);
    const std::string primary_dir = bam_info_dir + "/primary_metadata";
    const std::string primary_tar = bam_info_dir + "/primary_metadata.tar";
    const std::string primary_archive = primary_tar + ".zst";
    const std::string error_metadata = bam_info_dir + "/error.baminfo";
    const std::string error_archive = error_metadata + ".zst";
    bool ok = true;

    fprintf(stderr, "[BAM-INFO-TAR] archiving %s -> %s\n",
            primary_dir.c_str(), primary_tar.c_str());
    const std::string tar_command =
        "tar -C " + shell_quote(bam_info_dir) +
        " -cf " + shell_quote(primary_tar) + " primary_metadata";
    if (!execute_system_command(tar_command, "archive primary BAM metadata failed")) {
        ok = false;
    } else {
        fprintf(stderr, "[BAM-INFO-ZSTD] compressing %s -> %s\n",
                primary_tar.c_str(), primary_archive.c_str());
        const std::string command =
            "zstd -q -T" + std::to_string(zstd_threads) +
            " -19 -f --rm " + shell_quote(primary_tar) +
            " -o " + shell_quote(primary_archive);
        if (!execute_system_command(command, "compress primary metadata tar with zstd -19 failed")) {
            ok = false;
        } else {
            std::error_code remove_ec;
            fs::remove_all(primary_dir, remove_ec);
            if (remove_ec) {
                fprintf(stderr, "[BAM-INFO] failed to remove %s: %s\n",
                        primary_dir.c_str(), remove_ec.message().c_str());
                ok = false;
            }
        }
    }

    if (file_exists(error_metadata)) {
        fprintf(stderr, "[BAM-INFO-ZSTD] compressing %s -> %s\n",
                error_metadata.c_str(), error_archive.c_str());
        const std::string command =
            "zstd -q -T" + std::to_string(zstd_threads) +
            " -19 -f --rm " + shell_quote(error_metadata) +
            " -o " + shell_quote(error_archive);
        if (!execute_system_command(command, "compress error BAM metadata with zstd -19 failed")) {
            ok = false;
        }
    }

    std::ofstream manifest(bam_info_dir + "/bam_info_manifest.txt",
                           std::ios::out | std::ios::trunc);
    manifest
        << "version=2\n"
        << "header=header.sam\n"
        << "secondary_supplementary=discarded\n"
        << "unpaired_primary=unpaired_primary.bam\n"
        << "primary_metadata=primary_metadata.tar.zst\n"
        << "primary_metadata_layout=zstd-compressed tar containing primary_metadata/metadata_<partition>.bin\n"
        << "primary_metadata_alignment=sorted blocks; one metadata record per restored read\n"
        << "primary_record_order=partition/window sorted; original global BAM order is not retained\n"
        << "error_metadata=error.baminfo.zst\n"
        << "error_metadata_alignment=record-for-record with sibling error.fastq\n"
        << "metadata_fields=bam1_core,quality_missing,cigar,raw_aux\n"
        << "sequence_source=xzip\n"
        << "quality_source=xzip\n"
        << "qname_source=archive qnames directory\n"
        << "full_flags_source=primary metadata bam1_core.flag\n"
        << "binary_codec=zstd-19\n"
        << "bam_sam_codec=none\n"
        << "raw_binary_intermediates_retained=0\n";
    return ok && manifest.good();
}

enum class QualityCodec {
    ZstdLossless,
    WebpLossless4,
    WebpLossy
};

const char* quality_codec_name(QualityCodec codec) {
    switch (codec) {
        case QualityCodec::ZstdLossless: return "zstd-lossless";
        case QualityCodec::WebpLossless4: return "webp-lossless-4";
        case QualityCodec::WebpLossy: return "webp-lossy";
    }
    return "zstd-lossless";
}

bool write_quality_codec(const std::string& parent_dir, QualityCodec codec) {
    std::ofstream out(parent_dir + "/quality_codec.txt", std::ios::out | std::ios::trunc);
    if (!out) return false;
    out << quality_codec_name(codec) << '\n';
    return out.good();
}

bool write_xzip_manifest(const std::string& parent_dir,
                         QualityCodec codec,
                         int read_length,
                         bool paired_end,
                         bool preserve_qname,
                         bool bam_info_full) {
    std::ofstream out(parent_dir + "/xzip_manifest.txt", std::ios::out | std::ios::trunc);
    if (!out) return false;
    out << "version=4\n"
        << "archive_layout=tar-then-zstd\n"
        << "quality_codec=" << quality_codec_name(codec) << '\n'
        << "read_length=" << read_length << '\n'
        << "paired_end=" << (paired_end ? 1 : 0) << '\n'
        << "qname_mode=" << (preserve_qname ? "preserved" : "generated") << '\n'
        << "bam_info=" << (bam_info_full ? "full" : "none") << '\n'
        << "can_restore_fastq=1\n"
        << "can_restore_bam=" << (bam_info_full ? 1 : 0) << '\n';
    return out.good();
}

bool read_xzip_manifest(const std::string& parent_dir,
                        OL_PAR* options) {
    std::ifstream in(parent_dir + "/xzip_manifest.txt");
    if (!in) return false;
    std::string line;
    options->preserve_qname = 0;
    options->bam_info_full = 0;
    while (std::getline(in, line)) {
        const size_t equal = line.find('=');
        if (equal == std::string::npos) continue;
        const std::string key = line.substr(0, equal);
        const std::string value = line.substr(equal + 1);
        try {
            if (key == "read_length") options->read_length = std::stoi(value);
            else if (key == "paired_end") options->paired_end = std::stoi(value);
            else if (key == "qname_mode") options->preserve_qname = (value == "preserved") ? 1 : 0;
            else if (key == "bam_info") options->bam_info_full = (value == "full") ? 1 : 0;
        } catch (...) {
            return false;
        }
    }
    return options->read_length > 0 &&
           (options->paired_end == 0 || options->paired_end == 1);
}

bool decompress_legacy_zstd_files(const std::string& root_dir) {
    // Compatibility path for version-2 archives that stored one .zst per file.
    std::vector<fs::path> files;
    std::error_code ec;
    for (fs::recursive_directory_iterator it(root_dir, ec), end; it != end && !ec; it.increment(ec)) {
        if (!it->is_regular_file(ec)) continue;
        const fs::path& path = it->path();
        if (path.extension() == ".zst") files.push_back(path);
    }
    if (ec) {
        std::cerr << "[XZIP] cannot enumerate staging files: " << ec.message() << std::endl;
        return false;
    }

    for (size_t index = 0; index < files.size(); ++index) {
        const fs::path& file = files[index];
        std::error_code size_ec;
        const uint64_t measured_input_bytes = fs::file_size(file, size_ec);
        const uint64_t input_bytes = size_ec ? 0 : measured_input_bytes;
        const auto start = PerfClock::now();
        fprintf(stderr,
                "[XZIP-ZSTD] %s %zu/%zu start file=%s input_bytes=%" PRIu64 "\n",
                "decompress",
                index + 1, files.size(), file.string().c_str(),
                input_bytes);

        const std::string command =
            "zstd -q -d -f --rm " + shell_quote(file.string());
        if (system(command.c_str()) != 0) {
            std::cerr << "[XZIP] zstd failed for " << file << std::endl;
            return false;
        }

        const fs::path output =
            fs::path(file.string().substr(0, file.string().size() - 4));
        size_ec.clear();
        const uint64_t output_bytes = fs::file_size(output, size_ec);
        const double seconds = elapsed_ms(start, PerfClock::now()) / 1000.0;
        fprintf(stderr,
                "[XZIP-ZSTD] %s %zu/%zu done file=%s output_bytes=%" PRIu64
                " ratio=%.4f elapsed=%.2f s\n",
                "decompress",
                index + 1, files.size(), output.string().c_str(),
                size_ec ? 0 : output_bytes,
                input_bytes > 0 ? static_cast<double>(output_bytes) / input_bytes : 0.0,
                seconds);
    }

    std::cout << "[XZIP] decompressed " << files.size()
              << " legacy per-file zstd entries" << std::endl;
    return true;
}

bool create_xzip_archive(const std::string& parent_dir,
                         const std::string& archive_path,
                         QualityCodec codec,
                         int thread_n,
                         bool preserve_qname,
                         bool bam_info_full) {
    // Clean up the staging directory used by the previous archive implementation.
    std::error_code legacy_stage_ec;
    fs::remove_all(archive_path + ".stage.tmp", legacy_stage_ec);

    std::vector<std::string> entries = {
        "xzip_manifest.txt", "quality_codec.txt", "pre_window_id.txt",
        "partition_read_offsets.txt", "final_store_with_name", "diff_base",
        "diff_seq", "byte_flags"
    };
    entries.push_back(codec == QualityCodec::ZstdLossless
        ? "quality_score" : "store_2_zero_lossless");
    if (preserve_qname) entries.push_back("qnames");
    if (bam_info_full) entries.push_back("bam_info");

    for (const std::string& entry : entries) {
        std::error_code ec;
        if (!fs::exists(fs::path(parent_dir) / entry, ec)) {
            std::cerr << "[XZIP] required entry is missing: " << entry << std::endl;
            return false;
        }
    }
    if (file_exists(parent_dir + "/error.fastq")) entries.push_back("error.fastq");

    const std::string file_list_path = archive_path + ".files.txt";
    const std::string tar_path = archive_path + ".tar";
    const std::string archive_tmp_path = archive_path + ".tmp";
    {
        std::ofstream file_list(file_list_path, std::ios::out | std::ios::trunc);
        if (!file_list) {
            std::cerr << "[XZIP] cannot create tar file list: " << file_list_path << std::endl;
            return false;
        }
        for (const std::string& entry : entries) file_list << entry << '\n';
    }

    std::error_code remove_ec;
    fs::remove(tar_path, remove_ec);
    remove_ec.clear();
    fs::remove(archive_tmp_path, remove_ec);

    std::cout << "[XZIP-TAR] archiving " << entries.size()
              << " required entries into " << tar_path << std::endl;
    const std::string tar_command =
        "tar -C " + shell_quote(parent_dir) +
        " -cf " + shell_quote(tar_path) +
        " -T " + shell_quote(file_list_path);
    if (!execute_system_command(tar_command, "create temporary tar archive failed")) {
        return false;
    }

    std::error_code size_ec;
    const uint64_t tar_bytes = fs::file_size(tar_path, size_ec);
    const int zstd_threads = std::max(1, thread_n);
    fprintf(stderr,
            "[XZIP-ZSTD] compress tar start file=%s input_bytes=%" PRIu64
            " threads=%d level=19\n",
            tar_path.c_str(), size_ec ? 0 : tar_bytes, zstd_threads);
    const auto zstd_start = PerfClock::now();
    const std::string zstd_command =
        "zstd -q -T" + std::to_string(zstd_threads) +
        " -19 -f " + shell_quote(tar_path) +
        " -o " + shell_quote(archive_tmp_path);
    if (!execute_system_command(zstd_command, "compress tar archive with zstd -19 failed")) {
        std::cerr << "[XZIP] retained partial archive for debugging: "
                  << archive_tmp_path << std::endl;
        return false;
    }

    size_ec.clear();
    const uint64_t archive_bytes = fs::file_size(archive_tmp_path, size_ec);
    if (size_ec) {
        std::cerr << "[XZIP] cannot inspect compressed archive: "
                  << size_ec.message() << std::endl;
        return false;
    }
    remove_ec.clear();
    fs::rename(archive_tmp_path, archive_path, remove_ec);
    if (remove_ec) {
        std::cerr << "[XZIP] cannot publish final archive: "
                  << remove_ec.message() << std::endl;
        return false;
    }
    fprintf(stderr,
            "[XZIP-ZSTD] compress tar done file=%s output_bytes=%" PRIu64
            " ratio=%.4f elapsed=%.2f s\n",
            archive_path.c_str(), size_ec ? 0 : archive_bytes,
            tar_bytes > 0 ? static_cast<double>(archive_bytes) / tar_bytes : 0.0,
            elapsed_ms(zstd_start, PerfClock::now()) / 1000.0);
    fs::remove(tar_path, remove_ec);
    remove_ec.clear();
    fs::remove(file_list_path, remove_ec);
    return true;
}

bool extract_xzip_archive(const std::string& archive_path,
                          const std::string& output_dir) {
    if (!file_exists(archive_path)) {
        std::cerr << "[XZIP] archive does not exist: " << archive_path << std::endl;
        return false;
    }
    std::error_code ec;
    fs::remove_all(output_dir, ec);
    ec.clear();
    if (!fs::create_directories(output_dir, ec) && ec) {
        std::cerr << "[XZIP] cannot create extraction directory: " << output_dir << std::endl;
        return false;
    }
    if (is_zstd_frame(archive_path)) {
        std::cout << "[XZIP] detected tar-then-zstd archive" << std::endl;
        const std::string command =
            "zstd -d -c " + shell_quote(archive_path) +
            " | tar -C " + shell_quote(output_dir) + " -xf -";
        return execute_system_command(command, "extract tar-then-zstd .xzip archive failed");
    }

    std::cout << "[XZIP] detected legacy per-file-zstd-in-tar archive" << std::endl;
    const std::string command =
        "tar -C " + shell_quote(output_dir) + " -xf " + shell_quote(archive_path);
    if (!execute_system_command(command, "extract legacy .xzip tar archive failed")) return false;
    return decompress_legacy_zstd_files(output_dir);
}

uint64_t peak_memory_bytes() {
    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#ifdef __APPLE__
    return static_cast<uint64_t>(usage.ru_maxrss);
#else
    return static_cast<uint64_t>(usage.ru_maxrss) * 1024;
#endif
}

bool write_partition_read_offsets(const std::string& parent_dir,
                                  const std::vector<uint64_t>& partition_block_counts) {
    std::ofstream out(parent_dir + "/partition_read_offsets.txt",
                      std::ios::out | std::ios::trunc);
    if (!out) return false;
    uint64_t offset = 0;
    for (uint64_t count : partition_block_counts) {
        out << offset << '\n';
        offset += count;
    }
    out << offset << '\n';
    return out.good();
}

bool cleanup_compression_intermediates(const std::string& work_dir) {
    const std::vector<std::string> generated_entries = {
        "tmp", "final_store_with_name", "diff_base", "diff_seq", "byte_flags",
        "quality_score", "store_2_zero_lossless", "quality_diff_base",
        "reverse_seq", "output_dir", "qnames", "bam_info", "quality_codec.txt",
        "pre_window_id.txt", "partition_read_offsets.txt", "xzip_manifest.txt",
        "error.fastq"
    };
    bool ok = true;
    for (const std::string& entry : generated_entries) {
        std::error_code ec;
        fs::remove_all(fs::path(work_dir) / entry, ec);
        if (ec) {
            ok = false;
            std::cerr << "[XZIP] failed to remove intermediate " << entry
                      << ": " << ec.message() << std::endl;
        }
    }
    return ok;
}

QualityCodec read_quality_codec(const std::string& parent_dir, bool cli_lossless) {
    std::ifstream in(parent_dir + "/quality_codec.txt");
    std::string codec;
    if (in >> codec) {
        if (codec == "webp-lossless-4") return QualityCodec::WebpLossless4;
        if (codec == "webp-lossy") return QualityCodec::WebpLossy;
        if (codec == "zstd-lossless") return QualityCodec::ZstdLossless;
        std::cerr << "[Quality] unknown codec marker: " << codec << std::endl;
    }
    return cli_lossless ? QualityCodec::ZstdLossless : QualityCodec::WebpLossy;
}

size_t count_observed_quality_scores() {
    size_t distinct = 0;
    for (const auto& seen : g_quality_score_seen) {
        if (seen.load(std::memory_order_relaxed)) ++distinct;
    }
    return distinct;
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
 void store_the_quality_score(const std::vector<std::string>& quality_score_blocks, int id, char *quality_score_dir)
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
    for (const std::string& str : quality_score_blocks)
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

    for (const auto& block : data.Final_blocks) {
        // 原有二进制写入逻辑
        std::vector<unsigned char> binary_data = convertStringToBits(block.huffman_code);
        uint8_t code_len = block.huffman_code.size();
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

    // FASTQ reconstruction only needs the reverse-strand bit.
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

    // qnames remain outside the .xzip archive for optional BAM restoration.
    std::string parent_dir = get_parent_dir(pos_dir);
    std::string flags_dir = parent_dir + "/qnames";
    if (!create_directory(flags_dir)) {
        return;  // 文件夹创建失败，直接返回
    }
    // Each partition stores one qname text file.
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

 void count_ACGT(const std::vector<Compress_bio_string_block>& bio_string_blocks, int id, char *pos_dir)
 {
     Compress_final_block_with_huffman_table final_block_with_huffman_table;
         std::unordered_map<char,uint64_t>freqMap;
         for (const Compress_bio_string_block& block : bio_string_blocks)
         {
             for (char c : block.bio_string) {
             freqMap[c]++;
         }
     }
     HuffmanCoder coder;
     auto root = coder.buildHuffmanTree(freqMap);
     coder.buildCodes(root, "", final_block_with_huffman_table.huffmanCode);
     final_block_with_huffman_table.Final_blocks.reserve(bio_string_blocks.size());
     for (const Compress_bio_string_block& block : bio_string_blocks)
     {
         Compress_final_block final_block;
         std::string encoded = coder.encode(block.bio_string, final_block_with_huffman_table.huffmanCode);
         final_block.huffman_code = std::move(encoded);
         int length_qname = block.lqname;
         memcpy(final_block.qname, block.qname, length_qname);
         final_block.lqname = length_qname;
         final_block.real_seq1 = block.real_seq1;
         final_block.real_seq2 = block.real_seq2;
         final_block.hn = block.hn;
         final_block.hn2 = block.hn2;
         final_block_with_huffman_table.Final_blocks.push_back(std::move(final_block));
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
        char *pos_dir, char *quality_score_dir,
        std::vector<Compress_block> blocks) {

    if (blocks.empty()) return;

    std::vector<std::string>                    quality_score_blocks;
    std::vector<Compress_bio_string_block>      bio_string_blocks;
    std::vector<std::pair<uint16_t, uint16_t>>  flag_blocks;
    std::vector<std::string>                    qname_blocks;
    std::vector<std::string>                    diff_seq_blocks;
    std::vector<std::string>                    diff_base_blocks;
    const size_t estimated_reads = blocks.size() * 2;
    quality_score_blocks.reserve(estimated_reads);
    bio_string_blocks.reserve(blocks.size());
    flag_blocks.reserve(blocks.size());
    qname_blocks.reserve(blocks.size());
    diff_seq_blocks.reserve(estimated_reads);
    diff_base_blocks.reserve(estimated_reads);

    const std::string parent_dir = get_parent_dir(pos_dir);
    if (g_global_o && g_global_o->bam_info_full) {
        const std::string metadata_dir = parent_dir + "/bam_info/primary_metadata";
        std::error_code metadata_ec;
        fs::create_directories(metadata_dir, metadata_ec);
        const std::string metadata_path =
            metadata_dir + "/metadata_" + std::to_string(id) + ".bin";
        if (metadata_ec || !write_bam_metadata_file(metadata_path, blocks)) {
            fprintf(stderr, "[BAM-INFO] failed to write primary metadata: %s\n",
                    metadata_path.c_str());
        } else {
            fprintf(stderr, "[BAM-INFO] wrote primary metadata: %s (%zu blocks)\n",
                    metadata_path.c_str(), blocks.size());
        }
    }

    for (Compress_block& block : blocks) {
        const bool block_is_paired = block.is_paired != 0;

        quality_score_blocks.push_back(std::move(block.quality_score1));
        if (block_is_paired) quality_score_blocks.push_back(std::move(block.quality_score2));
        flag_blocks.push_back(std::make_pair(block.flags_1, block.flags_2));
        qname_blocks.emplace_back(block.qname, block.lqname);
        diff_seq_blocks.push_back(std::move(block.diff_seq1));
        if (block_is_paired) diff_seq_blocks.push_back(std::move(block.diff_seq2));
        diff_base_blocks.push_back(std::move(block.diff_base1));
        if (block_is_paired) diff_base_blocks.push_back(std::move(block.diff_base2));

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
        int length_qname = block.lqname;
        memcpy(bio_string_block.qname, block.qname, length_qname);
        bio_string_block.lqname    = length_qname;
        bio_string_block.real_seq1 = std::move(block.real_seq1);
        bio_string_block.real_seq2 = std::move(block.real_seq2);
        bio_string_block.hn2 = block_is_paired ? 1 : 0;
        bio_string_blocks.push_back(std::move(bio_string_block));

        prev_window_id = window_id1; // 更新，供下一条使用
    }

    store_the_quality_score(quality_score_blocks, id, quality_score_dir);
    count_ACGT(bio_string_blocks, id, pos_dir);
    store_the_flag(flag_blocks, pos_dir, id);
    if (g_global_o && g_global_o->preserve_qname) {
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
        if (!in_file.read(reinterpret_cast<char*>(&code_len), sizeof(uint8_t)))
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
        if (!in_file.read(reinterpret_cast<char*>(&code_len), sizeof(uint16_t)))
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

    return true;
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
    if (original_bits_length == 0) {
        std::cout << "diff_seq_read: 差异碱基为空（完全匹配分区）。" << std::endl;
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
 uint64_t reverse_dev_num(std::string dev_string)
 {
     uint64_t value = 0;
     for (char base : dev_string) {
         uint8_t bits;
         switch (base) {
             case 'A': bits = 0; break;
             case 'C': bits = 1; break;
             case 'G': bits = 2; break;
             case 'T': bits = 3; break;
             default: continue;
         }
         value = (value << 2) | bits;
     }
     return value;
 }

std::vector<Compress_block> bio_string_to_num(
    std::vector<Compress_bio_string_block> bio_string_blocks,
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
        pre_window_id = window_id1; // 更新，供下一条使用
        blocks.push_back(std::move(reverse_block));
    }

    return blocks;
}
int create_output_dir_c(char* webp_dir, char* output_dir, int output_dir_len) {
    if (webp_dir == NULL || output_dir == NULL || output_dir_len <= 0) {
        fprintf(stderr, "invalid output directory arguments\n");
        return -1;
    }

    const std::string path =
        (fs::path(webp_dir).parent_path() / "output_dir").string();
    if (path.size() + 1 > static_cast<size_t>(output_dir_len)) {
        fprintf(stderr, "output directory path is too long\n");
        return -1;
    }

    std::error_code ec;
    fs::create_directories(path, ec);
    if (ec) {
        fprintf(stderr, "cannot create output directory %s: %s\n",
                path.c_str(), ec.message().c_str());
        return -1;
    }
    memcpy(output_dir, path.c_str(), path.size() + 1);
    return 0;
}

struct RawBAMData {
    bam1_t* bam_ptr;  // 原始BAM数据（消费者处理后归还对象池）
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

        void push_batch(std::vector<RawBAMData>&& batch) {
            std::unique_lock<std::mutex> lock(mutex_);
            prod_cond_.wait(lock, [this]() { return queue_.size() < max_size_ || is_stop_; });
            if (is_stop_) return;
            queue_.push(std::move(batch));
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

// 4. 全局共享变量（线程安全）
const size_t MAX_QUEUE_SIZE = 32;
SafeTaskQueue g_task_queue(MAX_QUEUE_SIZE);          // 任务队列（生产者→消费者）
std::atomic<uint64_t> g_success_reads{0};      // 成功进入主压缩 block 的 reads
std::atomic<uint64_t> g_error_reads{0};        // 写入 error.fastq 的 reads
std::atomic<uint64_t> g_producer_primary_reads{0};   // producer 读到的 primary reads
std::atomic<uint64_t> g_producer_paired_reads{0};    // producer 成功配对送入 queue 的 reads
std::atomic<uint64_t> g_producer_unpaired_reads{0};  // producer 没配上的 primary reads
std::atomic<uint64_t> g_skipped_secondary{0};
std::atomic<uint64_t> g_skipped_supplementary{0};
std::atomic<uint64_t> g_reject_unmapped{0};
std::atomic<uint64_t> g_reject_read_length{0};
std::atomic<uint64_t> g_reject_invalid_tid{0};
std::atomic<uint64_t> g_reject_invalid_position{0};
bool g_input_is_paired = true;
CLASSIFY_SHARE_DATA* g_share = nullptr;

OL_PAR* g_global_o = nullptr;         // 全局参数（只读，线程安全）
const size_t BATCH_SIZE = 32768;
const size_t WRITE_THRESHOLD = 32768;

class BamRecordPool {
private:
    static constexpr size_t REFILL_SIZE = 4096;
    const size_t max_cached_;
    std::mutex mutex_;
    std::vector<bam1_t*> shared_cache_;
    std::vector<bam1_t*> producer_cache_;
    std::atomic<uint64_t> created_{0};
    std::atomic<uint64_t> reused_{0};
    std::atomic<uint64_t> returned_{0};
    std::atomic<uint64_t> destroyed_{0};

public:
    explicit BamRecordPool(size_t max_cached) : max_cached_(max_cached) {
        shared_cache_.reserve(max_cached_);
        producer_cache_.reserve(REFILL_SIZE);
    }

    ~BamRecordPool() {
        for (bam1_t* record : producer_cache_) bam_destroy1(record);
        for (bam1_t* record : shared_cache_) bam_destroy1(record);
    }

    bam1_t* acquire() {
        if (producer_cache_.empty()) {
            std::lock_guard<std::mutex> lock(mutex_);
            const size_t refill_count = std::min(REFILL_SIZE, shared_cache_.size());
            for (size_t i = 0; i < refill_count; ++i) {
                producer_cache_.push_back(shared_cache_.back());
                shared_cache_.pop_back();
            }
        }

        bam1_t* destination = nullptr;
        if (!producer_cache_.empty()) {
            destination = producer_cache_.back();
            producer_cache_.pop_back();
            reused_.fetch_add(1, std::memory_order_relaxed);
        } else {
            destination = bam_init1();
            if (!destination) {
                fprintf(stderr, "Fatal error: bam_init1 failed in BAM record pool\n");
                exit(EXIT_FAILURE);
            }
            created_.fetch_add(1, std::memory_order_relaxed);
        }
        return destination;
    }

    void release_one(bam1_t* record) {
        if (!record) return;
        std::vector<bam1_t*> records{record};
        release_batch(records);
    }

    void release_batch(std::vector<bam1_t*>& records) {
        if (records.empty()) return;

        size_t cache_count = 0;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            cache_count = std::min(records.size(), max_cached_ - shared_cache_.size());
            shared_cache_.insert(shared_cache_.end(), records.begin(), records.begin() + cache_count);
        }

        returned_.fetch_add(cache_count, std::memory_order_relaxed);
        for (size_t i = cache_count; i < records.size(); ++i) {
            bam_destroy1(records[i]);
        }
        destroyed_.fetch_add(records.size() - cache_count, std::memory_order_relaxed);
        records.clear();
    }

    void print_stats() {
        std::lock_guard<std::mutex> lock(mutex_);
        fprintf(stderr,
                "[bam-pool] created=%" PRIu64 " reused=%" PRIu64
                " returned=%" PRIu64 " destroyed=%" PRIu64
                " cached=%zu producer_cached=%zu\n",
                created_.load(std::memory_order_relaxed),
                reused_.load(std::memory_order_relaxed),
                returned_.load(std::memory_order_relaxed),
                destroyed_.load(std::memory_order_relaxed),
                shared_cache_.size(),
                producer_cache_.size());
    }
};

BamRecordPool g_bam_record_pool(BATCH_SIZE * 2);
std::vector<std::unique_ptr<std::mutex>> g_partition_locks;
std::mutex g_fq_write_mutex;
std::string g_error_fastq_path;
FILE* g_error_fq_fp = nullptr;
FILE* g_error_baminfo_fp = nullptr;
std::atomic<uint64_t> g_error_baminfo_records{0};
samFile* g_unpaired_primary_bam = nullptr;
std::string g_bam_info_dir;
thread_local std::string tl_error_fq_buffer;
thread_local std::vector<char> tl_error_baminfo_buffer;
const size_t ERROR_FQ_FLUSH_THRESHOLD = 4 * 1024 * 1024;  // 4MB
static char g_error_fq_stdio_buffer[1 << 20];
std::atomic<uint64_t> g_error_fastq_placeholder_id{1};

void append_error_metadata(const bam1_t* read) {
    append_bam_record_metadata(tl_error_baminfo_buffer, read);
}

void init_bam_info_outputs(const bam_hdr_t* header) {
    const std::string parent_dir = get_parent_dir(g_share->o->pos_dir);
    g_bam_info_dir = parent_dir + "/bam_info";
    std::error_code ec;
    fs::remove_all(g_bam_info_dir + "/primary_metadata", ec);
    ec.clear();
    fs::remove(g_bam_info_dir + "/primary_metadata.tar", ec);
    ec.clear();
    fs::remove(g_bam_info_dir + "/primary_metadata.tar.zst", ec);
    ec.clear();
    fs::remove(g_bam_info_dir + "/error.baminfo.zst", ec);
    ec.clear();
    fs::create_directories(g_bam_info_dir + "/primary_metadata", ec);
    if (ec) {
        fprintf(stderr, "Error: cannot create bam_info directories: %s\n",
                ec.message().c_str());
        exit(EXIT_FAILURE);
    }

    const std::string header_path = g_bam_info_dir + "/header.sam";
    std::ofstream header_out(header_path, std::ios::binary | std::ios::trunc);
    if (header && header->text && header->l_text > 0) {
        header_out.write(header->text, header->l_text);
    }
    header_out.close();

    fs::remove(g_bam_info_dir + "/secondary_supplementary.bam", ec);
    ec.clear();
    const std::string unpaired_path = g_bam_info_dir + "/unpaired_primary.bam";
    g_unpaired_primary_bam = sam_open(unpaired_path.c_str(), "wb");
    if (!g_unpaired_primary_bam ||
        sam_hdr_write(g_unpaired_primary_bam, header) < 0) {
        fprintf(stderr, "Error: failed to initialize sibling BAM-info outputs\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "[BAM-INFO] initialized %s\n", g_bam_info_dir.c_str());
}

void close_special_bam_outputs() {
    if (g_unpaired_primary_bam) {
        sam_close(g_unpaired_primary_bam);
        g_unpaired_primary_bam = nullptr;
    }
}

void init_error_fastq_writer() {
    g_error_fastq_path = get_parent_dir(g_share->o->pos_dir) + "/error.fastq";
    g_error_fq_fp = fopen(g_error_fastq_path.c_str(), "w");
    if (!g_error_fq_fp) {
        fprintf(stderr, "Error: Failed to open error fq file %s\n", g_error_fastq_path.c_str());
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
    setvbuf(g_error_fq_fp, g_error_fq_stdio_buffer, _IOFBF, sizeof(g_error_fq_stdio_buffer));

    if (g_global_o && g_global_o->bam_info_full) {
        const std::string error_baminfo_path = g_bam_info_dir + "/error.baminfo";
        g_error_baminfo_fp = fopen(error_baminfo_path.c_str(), "wb");
        if (!g_error_baminfo_fp) {
            fprintf(stderr, "Error: Failed to open %s\n", error_baminfo_path.c_str());
            exit(EXIT_FAILURE);
        }
        static const char magic[8] = {'X', 'E', 'R', 'R', 'M', 'E', 'T', '1'};
        const uint32_t version = 2;
        const uint64_t record_count_placeholder = 0;
        fwrite(magic, sizeof(magic), 1, g_error_baminfo_fp);
        fwrite(&version, sizeof(version), 1, g_error_baminfo_fp);
        fwrite(&record_count_placeholder, sizeof(record_count_placeholder), 1, g_error_baminfo_fp);
    }
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
    const size_t metadata_written = g_error_baminfo_fp
        ? fwrite(tl_error_baminfo_buffer.data(), 1, tl_error_baminfo_buffer.size(),
                 g_error_baminfo_fp)
        : tl_error_baminfo_buffer.size();
    if (g_error_baminfo_fp && metadata_written != tl_error_baminfo_buffer.size()) {
        fprintf(stderr,
                "Error: Failed to fully write error BAM-info buffer, expected=%zu actual=%zu\n",
                tl_error_baminfo_buffer.size(), metadata_written);
    }

    if (force) {
        fflush(g_error_fq_fp);
        if (g_error_baminfo_fp) fflush(g_error_baminfo_fp);
    }
    tl_error_fq_buffer.clear();
    tl_error_baminfo_buffer.clear();
}

void close_error_fastq_writer() {
    if (g_error_fq_fp) {
        fflush(g_error_fq_fp);
        fclose(g_error_fq_fp);
        g_error_fq_fp = nullptr;
    }
    if (g_error_baminfo_fp) {
        const uint64_t record_count =
            g_error_baminfo_records.load(std::memory_order_relaxed);
        fflush(g_error_baminfo_fp);
        fseek(g_error_baminfo_fp, sizeof(char) * 8 + sizeof(uint32_t), SEEK_SET);
        fwrite(&record_count, sizeof(record_count), 1, g_error_baminfo_fp);
        fclose(g_error_baminfo_fp);
        g_error_baminfo_fp = nullptr;
    }
    if (g_global_o && g_global_o->bam_info_full) {
        fprintf(stderr, "[BAM-INFO] error.baminfo records=%" PRIu64 "\n",
                g_error_baminfo_records.load(std::memory_order_relaxed));
    }
}

template <typename T>
void append_binary(std::vector<char>& buffer, const T& value) {
    const char* bytes = reinterpret_cast<const char*>(&value);
    buffer.insert(buffer.end(), bytes, bytes + sizeof(T));
}

void append_binary_data(std::vector<char>& buffer, const void* data, size_t size) {
    if (size == 0) return;
    const char* bytes = reinterpret_cast<const char*>(data);
    buffer.insert(buffer.end(), bytes, bytes + size);
}

void append_binary_string(std::vector<char>& buffer, const std::string& value) {
    const size_t size = value.size();
    append_binary(buffer, size);
    append_binary_data(buffer, value.data(), size);
}

void serialize_block(std::vector<char>& buffer, const Compress_block& block) {
    append_binary(buffer, block.window_id);
    append_binary(buffer, block.hap_offset);
    append_binary_data(buffer, block.qname, sizeof(block.qname));
    append_binary(buffer, block.lqname);
    append_binary(buffer, block.flags_1);
    append_binary_string(buffer, block.quality_score1);
    append_binary_string(buffer, block.real_seq1);
    append_binary_string(buffer, block.diff_seq1);
    append_binary_string(buffer, block.diff_base1);
    append_bam_metadata(buffer, block.bam_meta1);

    append_binary(buffer, block.window_id2);
    append_binary(buffer, block.hap_offset2);
    append_binary_data(buffer, block.qname2, sizeof(block.qname2));
    append_binary(buffer, block.lqname2);
    append_binary(buffer, block.flags_2);
    append_binary(buffer, block.is_paired);
    append_binary_string(buffer, block.quality_score2);
    append_binary_string(buffer, block.real_seq2);
    append_binary_string(buffer, block.diff_seq2);
    append_binary_string(buffer, block.diff_base2);
    append_bam_metadata(buffer, block.bam_meta2);
    append_binary(buffer, block.error_flag);
}

void flush_local_blocks(const std::vector<Compress_block>& blocks) {
    std::unordered_map<int, std::vector<const Compress_block*>> by_partition;
    by_partition.reserve(std::min<size_t>(blocks.size(), g_share->blockFiles.size()));

    for (const auto& block : blocks) {
        const int partition_idx = block.window_id / 10000;
        if (partition_idx < 0 ||
            partition_idx >= static_cast<int>(g_share->blockFiles.size())) continue;
        by_partition[partition_idx].push_back(&block);
    }

    for (auto& [partition_idx, partition_blocks] : by_partition) {
        FILE* fp = g_share->blockFiles[partition_idx];
        if (!fp) continue;

        std::vector<char> output;
        output.reserve(partition_blocks.size() *
                       static_cast<size_t>(2 * g_global_read_length + 256));
        for (const Compress_block* block : partition_blocks) {
            serialize_block(output, *block);
        }

        std::lock_guard<std::mutex> lock(*g_partition_locks[partition_idx]);
        const size_t written = fwrite(output.data(), 1, output.size(), fp);
        if (written != output.size()) {
            fprintf(stderr,
                    "Error: partition %d short write, expected=%zu actual=%zu\n",
                    partition_idx, output.size(), written);
        }
    }
}

// 生产者线程：读取BAM文件，生成任务包
void producer_thread(htsFile* input_file, bam_hdr_t* header) {
    bam1_t* read = g_bam_record_pool.acquire();
    std::vector<RawBAMData> current_batch;
    current_batch.reserve(BATCH_SIZE);

    bam1_t* r1 = nullptr;
    ProducerPerf perf;
    bool is_first_in_batch = true;
    PerfClock::time_point batch_start;
    size_t unpaired_warnings = 0;

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
            current_batch.push_back({read});
            read = g_bam_record_pool.acquire();
            g_producer_paired_reads++;
            perf.read_pairs++;

            if (current_batch.size() >= BATCH_SIZE) {
                auto before_push = PerfClock::now();
                perf.batch_fill_ms += elapsed_ms(batch_start, before_push);

                auto push_start = PerfClock::now();
                g_task_queue.push_batch(std::move(current_batch));
                auto push_end = PerfClock::now();
                perf.push_ms += elapsed_ms(push_start, push_end);

                perf.batch_cnt++;
                current_batch = {};
                current_batch.reserve(BATCH_SIZE);
                is_first_in_batch = true;
            }
            continue;
        }

        if (r1 == nullptr) {
            r1 = read;
            read = g_bam_record_pool.acquire();
        } else {
            if (strcmp(read_name, bam_get_qname(r1)) == 0) {
                current_batch.push_back({r1});
                current_batch.push_back({read});
                perf.read_pairs++;
                g_producer_paired_reads += 2;

                r1 = nullptr;
                read = g_bam_record_pool.acquire();

                if (current_batch.size() >= BATCH_SIZE) {
                    auto before_push = PerfClock::now();
                    perf.batch_fill_ms += elapsed_ms(batch_start, before_push);

                    auto push_start = PerfClock::now();
                    g_task_queue.push_batch(std::move(current_batch));
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

                    current_batch = {};
                    current_batch.reserve(BATCH_SIZE);
                    is_first_in_batch = true;
                }
            } else {
                g_producer_unpaired_reads++;

                if (unpaired_warnings++ < 10) {
                    fprintf(stderr,
                            "Warning: unpaired/non-adjacent primary read saved separately: %s\n",
                            bam_get_qname(r1));
                }
                if (g_global_o && g_global_o->bam_info_full &&
                    sam_write1(g_unpaired_primary_bam, header, r1) < 0) {
                    fprintf(stderr, "Error: failed writing unpaired primary BAM record\n");
                }
                g_bam_record_pool.release_one(r1);
                r1 = read;
                read = g_bam_record_pool.acquire();
            }
        }
    }

    if (r1 != nullptr) {
        g_producer_unpaired_reads++;

        if (unpaired_warnings++ < 10) {
            fprintf(stderr, "Warning: Unpaired read saved separately: %s\n", bam_get_qname(r1));
        }
        if (g_global_o && g_global_o->bam_info_full &&
            sam_write1(g_unpaired_primary_bam, header, r1) < 0) {
            fprintf(stderr, "Error: failed writing final unpaired primary BAM record\n");
        }
        g_bam_record_pool.release_one(r1);
        r1 = nullptr;
    }

    if (!current_batch.empty()) {
        auto before_push = PerfClock::now();
        perf.batch_fill_ms += elapsed_ms(batch_start, before_push);

        auto push_start = PerfClock::now();
        g_task_queue.push_batch(std::move(current_batch));
        auto push_end = PerfClock::now();
        perf.push_ms += elapsed_ms(push_start, push_end);

        perf.batch_cnt++;
    }

    g_task_queue.stop();
    g_bam_record_pool.release_one(read);

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

void write_error_reads_to_fq(const std::string& qname1, const std::string& seq1, const std::string& qual1,
                             const std::string& qname2, const std::string& seq2, const std::string& qual2,
                             const bam1_t* read1, const bam1_t* read2) {
    std::string out_qname1 = qname1;
    std::string out_qname2 = qname2;
    if (!g_global_o || !g_global_o->preserve_qname) {
        const uint64_t id = g_error_fastq_placeholder_id.fetch_add(1, std::memory_order_relaxed);
        out_qname1 = "error_" + std::to_string(id);
        out_qname2 = out_qname1;
    }
    const std::string fastq =
        "@" + out_qname1 + " 1:N:0:\n" + seq1 + "\n+\n" + qual1 + "\n" +
        "@" + out_qname2 + " 2:N:0:\n" + seq2 + "\n+\n" + qual2 + "\n";
    if (tl_error_fq_buffer.capacity() == 0) {
        tl_error_fq_buffer.reserve(ERROR_FQ_FLUSH_THRESHOLD + 1024);
        tl_error_baminfo_buffer.reserve(ERROR_FQ_FLUSH_THRESHOLD);
    }
    tl_error_fq_buffer += fastq;
    if (g_global_o && g_global_o->bam_info_full) {
        append_error_metadata(read1);
        append_error_metadata(read2);
        g_error_baminfo_records.fetch_add(2, std::memory_order_relaxed);
    }
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
    if (read->core.l_qseq > 0 && read->core.l_qseq < 400) {
        get_bam_seq(0, read->core.l_qseq, seq, read);
        uint8_t* qual = bam_get_qual(read);
        qual_str.clear();
        qual_str.reserve(read->core.l_qseq);
        for (int j = 0; j < read->core.l_qseq; ++j) {
            qual_str.push_back((char)(qual[j] + 33));
        }
    }

    if ((read->core.flag & BAM_FUNMAP) || read->core.n_cigar == 0) {
        g_reject_unmapped.fetch_add(1, std::memory_order_relaxed);
        return false;
    }
    if (read->core.l_qseq != g_global_read_length ||
        g_global_read_length <= 0 ||
        g_global_read_length >= 400) {
        g_reject_read_length.fetch_add(1, std::memory_order_relaxed);
        return false;
    }

    int chr_id = read->core.tid;
    if (chr_id < 0 || chr_id >= static_cast<int>(genome.size())) {
        g_reject_invalid_tid.fetch_add(1, std::memory_order_relaxed);
        return false;
    }

    uint wb_ref_pos = 0;
    const std::string& full_chrom_seq = genome[chr_id].second;
    uint local_offset = read->core.pos;
    if (local_offset == 0) window_id = chr_bg_wb_ID[chr_id];
    else window_id = chr_bg_wb_ID[chr_id] + (local_offset - 1) / g_global_read_length;
    local_offset -= (window_id - chr_bg_wb_ID[chr_id]) * g_global_read_length;
    if (local_offset >= static_cast<uint32_t>(g_global_read_length)) {
        window_id++;
        local_offset -= g_global_read_length;
    }
    wb_ref_pos = local_offset;

    uint32_t st_pos = read->core.pos + 1 - wb_ref_pos;
    uint32_t start_idx = st_pos - 1;
    if (start_idx >= full_chrom_seq.size()) {
        g_reject_invalid_position.fetch_add(1, std::memory_order_relaxed);
        return false;
    }

    int head_clip = get_cigar_head_clip_length(read->core.n_cigar, bam_get_cigar(read));
    int eo = (int)wb_ref_pos - head_clip;
    if (eo < 0 && window_id > 0 && window_id > chr_bg_wb_ID[chr_id]) {
        window_id--;
        eo += g_global_read_length;

        uint32_t new_st_pos = (window_id - chr_bg_wb_ID[chr_id]) * g_global_read_length + 1;
        uint32_t new_start_idx = new_st_pos - 1;
        if (new_start_idx >= full_chrom_seq.size()) {
            g_reject_invalid_position.fetch_add(1, std::memory_order_relaxed);
            return false;
        }
        start_idx = new_start_idx;
    } else if (eo < 0) {
        eo = 0;
    }
    effective_offset = (uint16_t)eo;

    diff_seq.clear();
    diff_base.clear();
    diff_seq.reserve(g_global_read_length);
    diff_base.reserve(g_global_read_length / 8);
    for (int k = 0; k < g_global_read_length; ++k) {
        int ref_k = eo + k;
        const int64_t genome_pos = static_cast<int64_t>(start_idx) + ref_k;
        bool mismatch = (ref_k < 0 ||
                         genome_pos < 0 ||
                         genome_pos >= static_cast<int64_t>(full_chrom_seq.size()))
                        ? true
                        : (seq[k] != full_chrom_seq[static_cast<size_t>(genome_pos)]);
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
                            int read_no,
                            const bam1_t* read) {
    std::string out_qname = qname;
    if (!g_global_o || !g_global_o->preserve_qname) {
        const uint64_t id = g_error_fastq_placeholder_id.fetch_add(1, std::memory_order_relaxed);
        out_qname = "error_" + std::to_string(id);
    }
    const std::string fastq =
        "@" + out_qname + (read_no == 2 ? " 2:N:0:\n" : " 1:N:0:\n") +
        seq + "\n+\n" + qual + "\n";
    if (tl_error_fq_buffer.capacity() == 0) {
        tl_error_fq_buffer.reserve(ERROR_FQ_FLUSH_THRESHOLD + 1024);
        tl_error_baminfo_buffer.reserve(ERROR_FQ_FLUSH_THRESHOLD);
    }
    tl_error_fq_buffer += fastq;
    if (g_global_o && g_global_o->bam_info_full) {
        append_error_metadata(read);
        g_error_baminfo_records.fetch_add(1, std::memory_order_relaxed);
    }
    flush_error_fq_buffer_if_needed(false);
}

// 7. 匹配相同read名称的逻辑
void consumer_thread() {
    thread_local std::vector<Compress_block> local_blocks;
    thread_local size_t block_count = 0;
    thread_local ConsumerPerf perf;
    thread_local std::array<bool, 256> quality_seen{};

    local_blocks.reserve(WRITE_THRESHOLD);

    const size_t SUMMARY_EVERY_BATCHES = 32;

    std::thread::id current_thread_id = std::this_thread::get_id();
    const size_t tid_hash = std::hash<std::thread::id>{}(current_thread_id);
    {
        std::lock_guard<std::mutex> lock(g_perf_print_mutex);
        fprintf(stderr, "[consumer-start %zu] ready\n", tid_hash);
    }

    std::vector<RawBAMData> raw_batch;
    bool first_batch = true;
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
        if (first_batch) {
            std::lock_guard<std::mutex> lock(g_perf_print_mutex);
            fprintf(stderr, "[consumer-first-batch %zu] records=%zu pairs=%zu\n",
                    tid_hash,
                    raw_batch.size(),
                    raw_batch.size() / (g_input_is_paired ? 2 : 1));
            first_batch = false;
        }

        if (g_input_is_paired && raw_batch.size() % 2 != 0) {
            fprintf(stderr, "Fatal error: Batch size is odd (not paired)\n");
            fflush(stderr);
            exit(EXIT_FAILURE);
        }

        const size_t reads_per_record = g_input_is_paired ? 2 : 1;
        uint64_t batch_success_reads = 0;
        uint64_t batch_error_reads = 0;
        std::vector<bam1_t*> completed_bams;
        completed_bams.reserve(raw_batch.size());

        for (size_t i = 0; i < raw_batch.size(); i += reads_per_record) {
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
                    write_error_reads_to_fq(name1, seq1, qual_str1, name2, seq2, qual_str2,
                                            read1, read2);
                    batch_error_reads += 2;
                } else {
                    write_error_read_to_fq(name1, seq1, qual_str1, 1, read1);
                    batch_error_reads += 1;
                }
                auto t_err_1 = PerfClock::now();
                perf.error_io_ms += elapsed_ms(t_err_0, t_err_1);
                perf.error_pair_cnt++;

                auto t_clean_0 = PerfClock::now();
                completed_bams.push_back(read1);
                if (read2) completed_bams.push_back(read2);
                auto t_clean_1 = PerfClock::now();
                perf.cleanup_ms += elapsed_ms(t_clean_0, t_clean_1);
                continue;
            }

            auto t_block_0 = PerfClock::now();
            for (unsigned char value : qual_str1) quality_seen[value] = true;
            for (unsigned char value : qual_str2) quality_seen[value] = true;

            if (read1->core.l_qname > sizeof(block.qname) ||
                (read2 && read2->core.l_qname > sizeof(block.qname2))) {
                auto t_err_0 = PerfClock::now();
                if (read2) {
                    write_error_reads_to_fq(name1, seq1, qual_str1, name2, seq2, qual_str2,
                                            read1, read2);
                    batch_error_reads += 2;
                } else {
                    write_error_read_to_fq(name1, seq1, qual_str1, 1, read1);
                    batch_error_reads += 1;
                }
                auto t_err_1 = PerfClock::now();
                perf.error_io_ms += elapsed_ms(t_err_0, t_err_1);
                perf.error_pair_cnt++;
                completed_bams.push_back(read1);
                if (read2) completed_bams.push_back(read2);
                continue;
            }

            block.window_id = window_id1;
            block.hap_offset = effective_offset1;
            memcpy(block.qname, name1, read1->core.l_qname);
            block.lqname = read1->core.l_qname;
            block.flags_1 = read1->core.flag;
            block.quality_score1 = std::move(qual_str1);
            block.real_seq1.assign(seq1, static_cast<size_t>(g_global_read_length));
            block.diff_seq1 = std::move(diff_seq1);
            block.diff_base1 = std::move(diff_base1);
            if (g_global_o && g_global_o->bam_info_full) {
                block.bam_meta1 = capture_bam_metadata(read1);
            }

            if (read2) {
                block.window_id2 = window_id2;
                block.hap_offset2 = effective_offset2;
                memcpy(block.qname2, name2, read2->core.l_qname);
                block.lqname2 = read2->core.l_qname;
                block.flags_2 = read2->core.flag;
                block.quality_score2 = std::move(qual_str2);
                block.real_seq2.assign(seq2, static_cast<size_t>(g_global_read_length));
                block.diff_seq2 = std::move(diff_seq2);
                block.diff_base2 = std::move(diff_base2);
                if (g_global_o && g_global_o->bam_info_full) {
                    block.bam_meta2 = capture_bam_metadata(read2);
                }
            }
            local_blocks.push_back(std::move(block));

            auto t_block_1 = PerfClock::now();
            perf.block_ms += elapsed_ms(t_block_0, t_block_1);

            block_count++;
            batch_success_reads += reads_per_record;

            if (block_count >= WRITE_THRESHOLD) {
                auto t_flush_0 = PerfClock::now();
                flush_local_blocks(local_blocks);
                auto t_flush_1 = PerfClock::now();
                perf.flush_ms += elapsed_ms(t_flush_0, t_flush_1);
                perf.flush_cnt++;

                local_blocks.clear();
                block_count = 0;
            }

            auto t_clean_0 = PerfClock::now();
            completed_bams.push_back(read1);
            if (read2) completed_bams.push_back(read2);
            auto t_clean_1 = PerfClock::now();
            perf.cleanup_ms += elapsed_ms(t_clean_0, t_clean_1);
        }

        auto t_pool_0 = PerfClock::now();
        g_bam_record_pool.release_batch(completed_bams);
        auto t_pool_1 = PerfClock::now();
        perf.cleanup_ms += elapsed_ms(t_pool_0, t_pool_1);

        g_success_reads.fetch_add(batch_success_reads, std::memory_order_relaxed);
        g_error_reads.fetch_add(batch_error_reads, std::memory_order_relaxed);

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
                    "flush=%.4f ms/pair cleanup=%.4f ms/pair error_pairs=%zu flush_cnt=%zu total_reads=%llu\n",
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
                    (unsigned long long)(g_success_reads.load(std::memory_order_relaxed) +
                                         g_error_reads.load(std::memory_order_relaxed)));
        }

        raw_batch.clear();
    }

    flush_error_fq_buffer_if_needed(true);

    if (!local_blocks.empty()) {
        auto t_flush_0 = PerfClock::now();
        flush_local_blocks(local_blocks);
        auto t_flush_1 = PerfClock::now();
        perf.flush_ms += elapsed_ms(t_flush_0, t_flush_1);
        perf.flush_cnt++;

        local_blocks.clear();
    }

    for (size_t i = 0; i < quality_seen.size(); ++i) {
        if (quality_seen[i]) {
            g_quality_score_seen[i].store(true, std::memory_order_relaxed);
        }
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
                "flush=%.4f ms/pair cleanup=%.4f ms/pair error_pairs=%zu flush_cnt=%zu total_reads=%llu\n",
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
                (unsigned long long)(g_success_reads.load(std::memory_order_relaxed) +
                                     g_error_reads.load(std::memory_order_relaxed)));
    }
}

bool compare_blocks(const Compress_block& a, const Compress_block& b) {
    if (a.window_id != b.window_id)
        return a.window_id < b.window_id;
    return a.hap_offset < b.hap_offset;  // 同 window 内按 wb_ref_pos 排序
}

bool read_bam_metadata(FILE* fp, BamReadMetadata& metadata, size_t& total_read) {
    total_read += fread(&metadata.core, sizeof(metadata.core), 1, fp);
    total_read += fread(&metadata.quality_missing, sizeof(metadata.quality_missing), 1, fp);
    uint32_t cigar_count = 0;
    total_read += fread(&cigar_count, sizeof(cigar_count), 1, fp);
    if (cigar_count > 1000000) return false;
    metadata.cigar.resize(cigar_count);
    if (cigar_count > 0) {
        total_read += fread(metadata.cigar.data(), sizeof(uint32_t), cigar_count, fp);
    }
    uint32_t aux_size = 0;
    total_read += fread(&aux_size, sizeof(aux_size), 1, fp);
    if (aux_size > 256 * 1024 * 1024u) return false;
    metadata.aux.resize(aux_size);
    if (aux_size > 0) total_read += fread(metadata.aux.data(), 1, aux_size, fp);
    return !ferror(fp);
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
        if (!read_bam_metadata(fp, block.bam_meta1, total_read)) break;

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
        if (!read_bam_metadata(fp, block.bam_meta2, total_read)) break;
        total_read += fread(&block.error_flag, sizeof(block.error_flag), 1, fp);
        if (total_read == 0) {
            break;
        }
        blocks.push_back(std::move(block));
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

    const size_t flush_size = 16 * 1024 * 1024;
    std::vector<char> output;
    output.reserve(flush_size + static_cast<size_t>(2 * g_global_read_length + 256));
    for (const auto& block : blocks) {
        serialize_block(output, block);
        if (output.size() >= flush_size) {
            fwrite(output.data(), 1, output.size(), fp);
            output.clear();
        }
    }
    if (!output.empty()) {
        fwrite(output.data(), 1, output.size(), fp);
    }

    fclose(fp);
}

void parallel_sort_and_write(
    std::vector<uint32_t>& pre_window_id_array,
    int total_groups,
    char* pos_dir,
    char* quality_dir,
    const std::string& work_dir,
    int max_threads,
    std::vector<uint64_t>* partition_block_counts = nullptr
) {
    const int worker_count = std::max(1, std::min(max_threads, total_groups));
    std::atomic<int> next_group{0};
    std::atomic<int> completed_groups{0};
    std::atomic<int> empty_groups{0};
    std::vector<std::thread> threads;
    threads.reserve(worker_count);

    for (int worker = 0; worker < worker_count; ++worker) {
        threads.emplace_back([&]() {
            while (true) {
                const int i = next_group.fetch_add(1, std::memory_order_relaxed);
                if (i >= total_groups) break;
            std::string file_name = work_dir + "/tmp/test." + std::to_string(i) + ".bin";

            std::vector<Compress_block> blocks = read_compress_blocks(file_name);
            if (blocks.empty()) {
                if (partition_block_counts && i < static_cast<int>(partition_block_counts->size())) {
                    (*partition_block_counts)[i] = 0;
                }
                fprintf(stderr, "Warning: No blocks in partition %d, skipping\n", i);
                empty_groups.fetch_add(1, std::memory_order_relaxed);
                completed_groups.fetch_add(1, std::memory_order_relaxed);
                continue;
            }

            uint32_t prev_window_id = pre_window_id_array[i];
            if (partition_block_counts && i < static_cast<int>(partition_block_counts->size())) {
                (*partition_block_counts)[i] = blocks.size();
            }
            sort_and_write_file(prev_window_id, i, pos_dir, quality_dir,
                                std::move(blocks));
            completed_groups.fetch_add(1, std::memory_order_relaxed);
            }
        });
    }

    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
    const int completed = completed_groups.load(std::memory_order_relaxed);
    const int empty = empty_groups.load(std::memory_order_relaxed);
    fprintf(stderr,
            "[partition-write-final] completed=%d total=%d empty=%d written=%d\n",
            completed, total_groups, empty, completed - empty);
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
    std::vector<uint32_t>& last_window_id
) {
    std::string file_name = work_dir + "/tmp/test." + std::to_string(partition_idx) + ".bin";
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

void generate_window_info(const std::vector<int>& chr_bg_wb_ID, int thread_num)
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

    // 3. 固定工作线程池处理染色体，避免复制整条参考序列和索引数组。
    int window_size = g_global_read_length;
    int read_length = 2 * g_global_read_length;

    int actual_thread_num = std::max(1, std::min(thread_num, total_chrs));
    std::cout << "实际使用线程数: " << actual_thread_num << "（总染色体数: " << total_chrs << "）" << std::endl;
    std::atomic<int> next_chr{0};
    std::vector<std::thread> threads;
    threads.reserve(actual_thread_num);
    for (int worker = 0; worker < actual_thread_num; ++worker) {
        threads.emplace_back([&] {
            while (true) {
                const int chr_id = next_chr.fetch_add(1, std::memory_order_relaxed);
                if (chr_id >= total_chrs) break;
                const auto& [start_window, end_window] = chr_window_ranges[chr_id];
                if (start_window > end_window) continue;
                process_window_range(chr_id, start_window, end_window,
                                     window_size, read_length,
                                     genome[chr_id].second, chr_bg_wb_ID);
            }
        });
    }
    for (auto& t : threads) {
        t.join();
    }

    std::cout << "所有染色体处理完成，总窗口数: " << g_window_info.size() << std::endl;
    return;
}

std::vector<int> build_chr_bg_wb_ID() {
    std::vector<int> chr_bg_wb_ID(genome.size() + 1, 0);
    for (size_t i = 0; i < genome.size(); ++i) {
        const std::string& current_seq = genome[i].second;
        int seq_len = current_seq.size();
        int window_cnt = (seq_len - g_global_read_length) / g_global_read_length + 1;
        chr_bg_wb_ID[i+1] = chr_bg_wb_ID[i] + window_cnt + 5;
    }
    return chr_bg_wb_ID;
}

int BWT_CLASSIFY_MAIN::init_run(int argc, char *argv[]){
    double cpu_time = cputime();
    double real_time1 = realtime();
    share = new CLASSIFY_SHARE_DATA();
    g_share = share;
    share->o = new OL_PAR();
    if (share->o->get_option(argc, argv) != 0) {
        delete share->o;
        delete share;
        return 2;
    }
    g_input_is_paired = (share->o->paired_end != 0);
    char *input_bam_fn = share->o->read_bam;
    htsFile *input_file = hts_open(input_bam_fn, "rb");
    if (!input_file) {
        fprintf(stderr, "Failed to open BAM file: %s\n", input_bam_fn);
        delete share->o;
        delete share;
        return 1;
    }
    const int total_thread_budget = std::max(1, share->o->thread_n);
    const int bam_io_threads = 1;
    bam_hdr_t *header = sam_hdr_read(input_file);
    share->header = header;
    int read_length = share->o->read_length;
    genome = read_fasta(std::string(share->o->fasta_path));
    g_global_o = share->o;
    g_global_read_length = read_length;

    chr_bg_wb_ID = build_chr_bg_wb_ID();
    uint64_t total_reference_bases = 0;
    for (const auto& chr : genome) total_reference_bases += chr.second.size();
    fprintf(stderr,
            "[Reference] sequences=%zu total_bases=%" PRIu64
            " read_length=%d windows=%d\n",
            genome.size(), total_reference_bases, g_global_read_length,
            chr_bg_wb_ID.empty() ? 0 : chr_bg_wb_ID.back());
    std::string work_dir = share->o->work_dir;
    std::error_code obsolete_flags_ec;
    fs::remove_all(fs::path(work_dir) / "flags", obsolete_flags_ec);
    if (obsolete_flags_ec) {
        fprintf(stderr, "[Flags] failed to remove obsolete full flags directory: %s\n",
                obsolete_flags_ec.message().c_str());
    }
    int window_file_partition = chr_bg_wb_ID[genome.size()] / 10000 + 1;  // 分区数量
    share->blockFiles.clear();  // 确保文件句柄列表为空
    for (int i = 0; i <= window_file_partition; ++i) {
        std::string file_name = work_dir + "/tmp/test." + std::to_string(i) + ".bin";
        FILE *fp = fopen(file_name.c_str(), "wb");  // 二进制写入模式
        if (!fp) {
            fprintf(stderr, "Failed to create partition file: %s\n", file_name.c_str());
            exit(EXIT_FAILURE);
        }
        share->blockFiles.emplace_back(fp);
    }

    g_partition_locks.clear();  // 确保为空
    for (int i = 0; i <= window_file_partition; ++i) {
        g_partition_locks.emplace_back(std::make_unique<std::mutex>());
    }
    fprintf(stderr, "Created %d partition files for window_id grouping\n", window_file_partition + 1);

    if (share->o->bam_info_full) {
        init_bam_info_outputs(header);
    }
    init_error_fastq_writer();

    // 4. 启动线程
    int thread_num = total_thread_budget;
    const int max_consumer_threads = 16;
    const int consumer_thread_num =
        std::max(1, std::min(max_consumer_threads, total_thread_budget - bam_io_threads));
    std::vector<std::thread> consumer_threads;

    generate_window_info(chr_bg_wb_ID, thread_num);

    // 启动消费者线程（N个）
    fprintf(stderr, "[Threads] total=%d bam_decode=%d consumers=%d\n",
            total_thread_budget, bam_io_threads, consumer_thread_num);
    fprintf(stderr,
            "[Producer] direct BAM ownership transfer enabled; "
            "HTSlib extra decode threads disabled\n");
    for (int i = 0; i < consumer_thread_num; i++) {
        consumer_threads.emplace_back(consumer_thread);
    }

    std::thread producer(
        producer_thread,       // 你的生产者线程函数
        input_file,            // BAM文件句柄（生产者负责读取）
        header                 // BAM头信息
    );

    producer.join();
    hts_close(input_file);
    if (share->o->bam_info_full) {
        close_special_bam_outputs();
    }

    // 9. 等待消费者线程结束（同原代码）
    for (auto& t : consumer_threads) {
        t.join();
    }
    g_bam_record_pool.print_stats();
    fprintf(stderr,
            "[Reject-final] unmapped_or_no_cigar=%" PRIu64
            " read_length_mismatch=%" PRIu64
            " invalid_tid=%" PRIu64
            " invalid_position=%" PRIu64
            " success_reads=%" PRIu64 " error_reads=%" PRIu64 "\n",
            g_reject_unmapped.load(std::memory_order_relaxed),
            g_reject_read_length.load(std::memory_order_relaxed),
            g_reject_invalid_tid.load(std::memory_order_relaxed),
            g_reject_invalid_position.load(std::memory_order_relaxed),
            g_success_reads.load(std::memory_order_relaxed),
            g_error_reads.load(std::memory_order_relaxed));

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
    const int partition_workers = std::max(1, std::min(max_threads, window_file_partition));
    std::atomic<int> next_partition{0};
    std::vector<std::thread> threads;
    threads.reserve(partition_workers);
    for (int worker = 0; worker < partition_workers; ++worker) {
        threads.emplace_back([&] {
            while (true) {
                const int i = next_partition.fetch_add(1, std::memory_order_relaxed);
                if (i >= window_file_partition) break;
                process_partition(i, work_dir, last_window_id);
            }
        });
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

    std::vector<uint64_t> partition_block_counts(window_file_partition + 1, 0);
    parallel_sort_and_write(
        pre_window_id,
        window_file_partition + 1,
        share->o->pos_dir,
        share->o->quality_score_dir,
        work_dir,
        share->o->thread_n,
        &partition_block_counts
    );
    if (share->o->bam_info_full) {
        if (!compress_bam_info_metadata(g_bam_info_dir, share->o->thread_n)) {
            std::cerr << "[BAM-INFO] metadata compression or manifest generation was incomplete"
                      << std::endl;
        } else {
            std::cout << "[BAM-INFO] BAM restoration data ready in "
                      << g_bam_info_dir << std::endl;
        }
    }

    std::string parent_dir_path = get_parent_dir(share->o->pos_dir);
    QualityCodec quality_codec = QualityCodec::WebpLossy;
    if (share->o->webp_lossless) {
        const size_t distinct_quality_scores = count_observed_quality_scores();
        if (distinct_quality_scores == 4) {
            quality_codec = QualityCodec::WebpLossless4;
            std::cout << "[Quality] detected " << distinct_quality_scores
                      << " quality values; encode with lossless WebP" << std::endl;
            webp_main_40(share->o->quality_score_dir, share->o->webp_dir, share->o->thread_n,
              chr_bg_wb_ID[genome.size()]/10000, g_global_read_length, 1, 100,
              share->o->webp_method, share->o->webp_width, share->o->webp_height);
        } else {
            quality_codec = QualityCodec::ZstdLossless;
            std::cout << "[Quality] detected " << distinct_quality_scores
                      << " quality values; package quality_score with zstd -19" << std::endl;
        }
    } else {
        quality_codec = QualityCodec::WebpLossy;
        std::cout << "[Quality] WebP lossy mode: encode quality scores with WebP" << std::endl;
        webp_main_40(share->o->quality_score_dir, share->o->webp_dir, share->o->thread_n,
          chr_bg_wb_ID[genome.size()]/10000, g_global_read_length, 0, share->o->webp_quality,
          share->o->webp_method, share->o->webp_width, share->o->webp_height);
    }
    if (!write_quality_codec(parent_dir_path, quality_codec)) {
        std::cerr << "[Quality] failed to write quality codec marker" << std::endl;
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
    if (!write_partition_read_offsets(parent_dir_path, partition_block_counts)) {
        std::cerr << "无法写入 partition_read_offsets.txt" << std::endl;
    }

    const std::string output_stem = fs::path(parent_dir_path).filename().string();
    const std::string archive_path = parent_dir_path + "/" + output_stem + ".xzip";
    const bool manifest_ok = write_xzip_manifest(
        parent_dir_path, quality_codec, g_global_read_length, g_input_is_paired,
        share->o->preserve_qname != 0, share->o->bam_info_full != 0);
    const bool archive_ok = manifest_ok && create_xzip_archive(
        parent_dir_path, archive_path, quality_codec, share->o->thread_n,
        share->o->preserve_qname != 0, share->o->bam_info_full != 0);
    if (archive_ok) {
        cleanup_compression_intermediates(parent_dir_path);
        std::error_code metrics_ec;
        fs::remove(parent_dir_path + "/" + output_stem + ".metrics.txt", metrics_ec);
        std::cout << "[XZIP] created " << archive_path << std::endl;
    } else {
        std::cerr << "[XZIP] packaging was incomplete; intermediate files were retained" << std::endl;
    }
    const double total_cpu_seconds = cputime() - cpu_time;
    const double total_real_seconds = realtime() - real_time1;

    // 8. 释放所有资源（同原代码）
    bam_hdr_destroy(header);
    delete share->o;
    delete share;



    fprintf(stderr, "Total CPU time: %.2f s\n", total_cpu_seconds);
    fprintf(stderr, "Total REAL time: %.2f s\n", total_real_seconds);
    fprintf(stderr, "Peak memory: %llu bytes\n",
            (unsigned long long)peak_memory_bytes());
    return archive_ok ? 0 : 1;
}


bool g_restore_bam_enabled = false;
std::string g_restore_bam_info_dir;
std::string g_restore_primary_metadata_dir;
std::string g_restore_error_metadata_path;
std::string g_restore_bam_output_path;

struct BamMetadataBlock {
    uint8_t is_paired = 0;
    BamReadMetadata read1;
    BamReadMetadata read2;
};

bool read_bam_metadata_stream(std::istream& in,
                              BamReadMetadata& metadata,
                              uint32_t format_version) {
    uint32_t cigar_count = 0;
    uint32_t aux_size = 0;
    if (!in.read(reinterpret_cast<char*>(&metadata.core), sizeof(metadata.core))) {
        return false;
    }
    if (format_version >= 2 &&
        !in.read(reinterpret_cast<char*>(&metadata.quality_missing),
                 sizeof(metadata.quality_missing))) {
        return false;
    }
    if (
        !in.read(reinterpret_cast<char*>(&cigar_count), sizeof(cigar_count)) ||
        cigar_count > 1000000) {
        return false;
    }
    metadata.cigar.resize(cigar_count);
    if (cigar_count > 0 &&
        !in.read(reinterpret_cast<char*>(metadata.cigar.data()),
                 cigar_count * sizeof(uint32_t))) {
        return false;
    }
    if (!in.read(reinterpret_cast<char*>(&aux_size), sizeof(aux_size)) ||
        aux_size > 256 * 1024 * 1024u) {
        return false;
    }
    metadata.aux.resize(aux_size);
    return aux_size == 0 ||
           static_cast<bool>(in.read(reinterpret_cast<char*>(metadata.aux.data()), aux_size));
}

class PrimaryMetadataReader {
public:
    explicit PrimaryMetadataReader(const std::string& path) : in_(path, std::ios::binary) {
        char magic[8] = {};
        uint32_t version = 0;
        if (!in_.read(magic, sizeof(magic)) ||
            std::memcmp(magic, "XBAMMET1", sizeof(magic)) != 0 ||
            !in_.read(reinterpret_cast<char*>(&version), sizeof(version)) ||
            (version != 1 && version != 2) ||
            !in_.read(reinterpret_cast<char*>(&block_count_), sizeof(block_count_))) {
            valid_ = false;
        } else {
            valid_ = true;
            version_ = version;
        }
    }

    bool valid() const { return valid_; }
    uint64_t block_count() const { return block_count_; }

    bool next(BamMetadataBlock& block) {
        if (!valid_ || records_read_ >= block_count_) return false;
        if (!in_.read(reinterpret_cast<char*>(&block.is_paired), sizeof(block.is_paired)) ||
            !read_bam_metadata_stream(in_, block.read1, version_) ||
            (block.is_paired && !read_bam_metadata_stream(in_, block.read2, version_))) {
            valid_ = false;
            return false;
        }
        ++records_read_;
        return true;
    }

private:
    std::ifstream in_;
    bool valid_ = false;
    uint64_t block_count_ = 0;
    uint64_t records_read_ = 0;
    uint32_t version_ = 0;
};

uint8_t bam_base_code(char base) {
    switch (base) {
        case '=': return 0;
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'M': case 'm': return 3;
        case 'G': case 'g': return 4;
        case 'R': case 'r': return 5;
        case 'S': case 's': return 6;
        case 'V': case 'v': return 7;
        case 'T': case 't': return 8;
        case 'W': case 'w': return 9;
        case 'Y': case 'y': return 10;
        case 'H': case 'h': return 11;
        case 'K': case 'k': return 12;
        case 'D': case 'd': return 13;
        case 'B': case 'b': return 14;
        default: return 15;
    }
}

bam1_t* rebuild_bam_record(const BamReadMetadata& metadata,
                           const std::string& qname,
                           const std::string& sequence,
                           const std::string& quality) {
    bam1_t* record = bam_init1();
    if (!record) return nullptr;

    record->core = metadata.core;
    record->core.n_cigar = static_cast<uint32_t>(metadata.cigar.size());
    record->core.l_qseq = static_cast<int32_t>(sequence.size());
    const size_t seq_bytes = (sequence.size() + 1) / 2;
    const size_t data_size = record->core.l_qname +
                             metadata.cigar.size() * sizeof(uint32_t) +
                             seq_bytes + sequence.size() + metadata.aux.size();
    record->data = static_cast<uint8_t*>(std::realloc(record->data, data_size));
    if (!record->data && data_size > 0) {
        bam_destroy1(record);
        return nullptr;
    }
    record->l_data = static_cast<int>(data_size);
    record->m_data = static_cast<uint32_t>(data_size);
    std::memset(record->data, 0, data_size);

    const size_t qname_capacity = record->core.l_qname > 0
        ? record->core.l_qname - 1 : 0;
    std::memcpy(record->data, qname.data(), std::min(qname.size(), qname_capacity));

    uint8_t* cigar_out = record->data + record->core.l_qname;
    if (!metadata.cigar.empty()) {
        std::memcpy(cigar_out, metadata.cigar.data(),
                    metadata.cigar.size() * sizeof(uint32_t));
    }
    uint8_t* seq_out = cigar_out + metadata.cigar.size() * sizeof(uint32_t);
    for (size_t i = 0; i < sequence.size(); ++i) {
        seq_out[i >> 1] |= bam_base_code(sequence[i]) << ((~i & 1) << 2);
    }
    uint8_t* qual_out = seq_out + seq_bytes;
    if (metadata.quality_missing) {
        std::memset(qual_out, 0xff, sequence.size());
    } else if (quality.size() == sequence.size()) {
        for (size_t i = 0; i < quality.size(); ++i) {
            qual_out[i] = static_cast<uint8_t>(
                std::max(0, static_cast<int>(static_cast<unsigned char>(quality[i])) - 33));
        }
    } else {
        std::memset(qual_out, 0xff, sequence.size());
    }
    if (!metadata.aux.empty()) {
        std::memcpy(qual_out + sequence.size(), metadata.aux.data(), metadata.aux.size());
    }
    return record;
}

std::string fastq_qname(const std::string& header) {
    size_t begin = (!header.empty() && header[0] == '@') ? 1 : 0;
    size_t end = header.find_first_of(" \t", begin);
    return header.substr(begin, end == std::string::npos ? end : end - begin);
}

bool prepare_bam_restore_inputs(const std::string& archive_path,
                                const std::string& unpacked_dir) {
    g_restore_bam_info_dir = unpacked_dir + "/bam_info";
    if (!directory_exists(g_restore_bam_info_dir) ||
        !file_exists(g_restore_bam_info_dir + "/header.sam") ||
        !file_exists(g_restore_bam_info_dir + "/primary_metadata.tar.zst")) {
        std::cerr << "[BAM-RESTORE] archive does not contain complete bam_info"
                  << std::endl;
        return false;
    }

    const std::string stage_dir = unpacked_dir + "/bam_restore_stage";
    std::error_code ec;
    fs::remove_all(stage_dir, ec);
    ec.clear();
    fs::create_directories(stage_dir, ec);
    if (ec) return false;

    const std::string extract_command =
        "zstd -q -d -c " + shell_quote(g_restore_bam_info_dir + "/primary_metadata.tar.zst") +
        " | tar -C " + shell_quote(stage_dir) + " -xf -";
    if (!execute_system_command(extract_command, "extract primary BAM metadata failed")) {
        return false;
    }
    g_restore_primary_metadata_dir = stage_dir + "/primary_metadata";

    g_restore_error_metadata_path.clear();
    const std::string error_archive = g_restore_bam_info_dir + "/error.baminfo.zst";
    if (file_exists(error_archive)) {
        g_restore_error_metadata_path = stage_dir + "/error.baminfo";
        const std::string error_command =
            "zstd -q -d -f " + shell_quote(error_archive) +
            " -o " + shell_quote(g_restore_error_metadata_path);
        if (!execute_system_command(error_command, "extract error BAM metadata failed")) {
            return false;
        }
    }

    if (g_global_o && g_global_o->restore_bam_path) {
        g_restore_bam_output_path = g_global_o->restore_bam_path;
    } else {
        const fs::path archive_parent = fs::absolute(fs::path(archive_path)).parent_path();
        const std::string stem = fs::path(archive_path).stem().string();
        g_restore_bam_output_path = (archive_parent / (stem + ".restored.bam")).string();
    }
    return true;
}

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

        next_blocks = bio_string_to_num(bio_string_blocks, prev_window_id);
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
        for(size_t j = 0; j < len1; ++j)
        {
            const unsigned char value = static_cast<unsigned char>(qual1[j]);
            if(value < 33) qual1[j] = '!';
            if(value > 126) qual1[j] = '~';
        }
        for(size_t j=0;j<qual2.size();j++)
        {
            const unsigned char value = static_cast<unsigned char>(qual2[j]);
            if(value < 33) qual2[j] = '!';
            if(value > 126) qual2[j] = '~';
        }
        quality_pairs.emplace_back(qual1, qual2);
    }
    return quality_pairs;
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

std::vector<uint64_t> read_partition_read_offsets(const std::string& base_dir) {
    std::vector<uint64_t> offsets;
    std::ifstream in(base_dir + "/partition_read_offsets.txt");
    std::string line;
    while (std::getline(in, line)) {
        try {
            offsets.push_back(static_cast<uint64_t>(std::stoull(line)));
        } catch (...) {}
    }
    return offsets;
}

std::vector<std::string> generate_qnames_for_partition(
    int id,
    size_t count,
    const std::vector<uint64_t>& partition_read_offsets) {
    std::vector<std::string> qnames;
    qnames.reserve(count);
    const uint64_t base = (id >= 0 && id < static_cast<int>(partition_read_offsets.size()))
        ? partition_read_offsets[id] : 0;
    for (size_t i = 0; i < count; ++i) {
        qnames.push_back(std::to_string(base + i + 1));
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
    std::ofstream& merged_r2,
    bool generate_qname,
    uint64_t& next_generated_id
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

        if (generate_qname) {
            const std::string name = std::to_string(next_generated_id++);
            h1 = "@" + name + "/1";
            h2 = "@" + name + "/2";
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
    std::ofstream& merged_r1,
    bool generate_qname,
    uint64_t& next_generated_id
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
        if (generate_qname) {
            h = "@" + std::to_string(next_generated_id++);
        }
        merged_r1 << h << "\n" << s << "\n" << p << "\n" << q << "\n";
        ++error_reads;
    }

    std::cout << "[MERGE] appended error.fastq reads = "
              << error_reads << std::endl;
}

bam_hdr_t* load_restore_bam_header() {
    std::ifstream in(g_restore_bam_info_dir + "/header.sam", std::ios::binary);
    if (!in.is_open()) return nullptr;
    const std::string text((std::istreambuf_iterator<char>(in)),
                           std::istreambuf_iterator<char>());
    return sam_hdr_parse(static_cast<int>(text.size()), text.c_str());
}

bool append_bam_file_records(const std::string& path,
                             samFile* output,
                             const bam_hdr_t* output_header,
                             uint64_t& written) {
    if (!file_exists(path)) return true;
    samFile* input = sam_open(path.c_str(), "rb");
    if (!input) return false;
    bam_hdr_t* input_header = sam_hdr_read(input);
    bam1_t* record = bam_init1();
    bool ok = input_header && record;
    while (ok && sam_read1(input, input_header, record) >= 0) {
        if (sam_write1(output, output_header, record) < 0) {
            ok = false;
            break;
        }
        ++written;
    }
    bam_destroy1(record);
    bam_hdr_destroy(input_header);
    sam_close(input);
    return ok;
}

bool append_error_bam_records(const std::string& error_fastq_path,
                              samFile* output,
                              const bam_hdr_t* output_header,
                              uint64_t& written) {
    if (g_restore_error_metadata_path.empty() ||
        !file_exists(g_restore_error_metadata_path) ||
        !file_exists(error_fastq_path)) {
        return true;
    }

    std::ifstream metadata(g_restore_error_metadata_path, std::ios::binary);
    std::ifstream fastq(error_fastq_path);
    char magic[8] = {};
    uint32_t version = 0;
    uint64_t record_count = 0;
    if (!metadata.read(magic, sizeof(magic)) ||
        std::memcmp(magic, "XERRMET1", sizeof(magic)) != 0 ||
        !metadata.read(reinterpret_cast<char*>(&version), sizeof(version)) ||
        (version != 1 && version != 2) ||
        !metadata.read(reinterpret_cast<char*>(&record_count), sizeof(record_count))) {
        return false;
    }

    for (uint64_t i = 0; i < record_count; ++i) {
        BamReadMetadata bam_metadata;
        std::string header, sequence, plus, quality;
        if (!read_bam_metadata_stream(metadata, bam_metadata, version) ||
            !std::getline(fastq, header) ||
            !std::getline(fastq, sequence) ||
            !std::getline(fastq, plus) ||
            !std::getline(fastq, quality)) {
            return false;
        }
        bam1_t* record = rebuild_bam_record(
            bam_metadata, fastq_qname(header), sequence, quality);
        if (!record || sam_write1(output, output_header, record) < 0) {
            bam_destroy1(record);
            return false;
        }
        bam_destroy1(record);
        ++written;
    }
    return true;
}

bool merge_restored_bam(int file_num,
                        const std::string& parent_dir,
                        const bam_hdr_t* header) {
    samFile* output = sam_open(g_restore_bam_output_path.c_str(), "wb");
    if (!output || sam_hdr_write(output, header) < 0) {
        if (output) sam_close(output);
        return false;
    }

    uint64_t primary_written = 0;
    uint64_t error_written = 0;
    uint64_t unpaired_written = 0;
    bool ok = true;
    const std::string parts_dir = parent_dir + "/restored_bam_parts";
    for (int id = 0; id <= file_num && ok; ++id) {
        ok = append_bam_file_records(
            parts_dir + "/part_" + std::to_string(id) + ".bam",
            output, header, primary_written);
    }
    if (ok) {
        ok = append_error_bam_records(parent_dir + "/error.fastq",
                                      output, header, error_written);
    }
    if (ok) {
        ok = append_bam_file_records(g_restore_bam_info_dir + "/unpaired_primary.bam",
                                     output, header, unpaired_written);
    }
    sam_close(output);

    if (ok) {
        std::error_code ec;
        fs::remove_all(parts_dir, ec);
        ec.clear();
        fs::remove_all(parent_dir + "/bam_restore_stage", ec);
        std::cout << "[BAM-RESTORE] wrote " << g_restore_bam_output_path
                  << " primary=" << primary_written
                  << " error=" << error_written
                  << " unpaired=" << unpaired_written << std::endl;
    }
    return ok;
}

void go_to_decompress_pos(std::vector<uint32_t> pre_window_id_vector,
    CLASSIFY_SHARE_DATA* share)
{
    double start_time = cputime();

    int read_length = share->o->read_length;

    constexpr bool SIMPLE_BLOCK_PROGRESS = true;
    constexpr int  BLOCK_LOG_INTERVAL = 10000;

    std::queue<Task> task_queue;
    std::mutex queue_mutex;
    std::condition_variable cv;
    std::atomic<bool> stop_flag(false);
    std::atomic<int> completed_tasks(0);
    std::atomic<bool> bam_restore_failed(false);

    int file_num = chr_bg_wb_ID[genome.size()] / 10000;
    const std::string decompress_parent_dir = get_parent_dir(share->o->pos_dir);
    const std::vector<uint64_t> partition_read_offsets =
        read_partition_read_offsets(decompress_parent_dir);
    if (!share->o->preserve_qname && partition_read_offsets.empty()) {
        std::cerr << "[FATAL] qname_mode=generated but partition_read_offsets.txt is missing"
                  << std::endl;
        return;
    }
    bam_hdr_t* restore_bam_header = nullptr;
    if (g_restore_bam_enabled) {
        restore_bam_header = load_restore_bam_header();
        if (!restore_bam_header) {
            std::cerr << "[BAM-RESTORE] cannot parse header.sam; disabling BAM restoration"
                      << std::endl;
            g_restore_bam_enabled = false;
        } else {
            std::error_code ec;
            fs::remove_all(get_parent_dir(share->o->pos_dir) + "/restored_bam_parts", ec);
            ec.clear();
            fs::create_directories(get_parent_dir(share->o->pos_dir) + "/restored_bam_parts", ec);
            if (ec) {
                std::cerr << "[BAM-RESTORE] cannot create part directory: "
                          << ec.message() << std::endl;
                g_restore_bam_enabled = false;
            }
        }
    }

    std::vector<int> task_ids_to_run;
    for (int i = 0; i <= file_num; ++i) {
        task_ids_to_run.push_back(i);
    }

    int total_tasks = static_cast<int>(task_ids_to_run.size());
    if (total_tasks == 0) {
        std::cerr << "[WARN] 没有可运行的 task，直接返回。" << std::endl;
        return;
    }

    std::cout << "[MAIN] file_num=" << file_num
              << ", total_tasks=" << total_tasks
              << ", read_length=" << read_length
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

            uint64_t required_diff_bases = 0;
            for (const std::string& diff_seq : diff_seqs) {
                required_diff_bases += static_cast<uint64_t>(
                    std::count(diff_seq.begin(), diff_seq.end(), '1'));
            }
            if (diff_base_decoded_str.size() != required_diff_bases) {
                err("diff_base 数量不匹配，decoded=" +
                    std::to_string(diff_base_decoded_str.size()) +
                    ", required=" + std::to_string(required_diff_bases),
                    task.file_idx);
                int done = ++completed_tasks;
                std::cout << "[TASK-END] task_idx=" << task.file_idx
                          << " status=diff_base_mismatch completed_tasks="
                          << done << "/" << total_tasks << std::endl;
                cv.notify_all();
                continue;
            }
            if (required_diff_bases == 0) {
                std::cout << "[LOAD] task_idx=" << task.file_idx
                          << " contains only exact-match reads; empty diff_base is valid"
                          << std::endl;
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
                    reconstruct_read(diff_seqs[process_idx], diff_base_decoded_str, diff_base_idx,
                                     block.window_id2, block.hap_offset2,
                                     reverse_read2, "R2", block_idx, fatal_error);
                    if (fatal_error) break;
                    process_idx++;
                }
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
            auto qnames = task.share->o->preserve_qname
                ? read_qnames(parent_dir, task.file_idx)
                : generate_qnames_for_partition(task.file_idx, reverse_pairs.size(),
                                                partition_read_offsets);
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

            if (g_restore_bam_enabled) {
                const std::string metadata_path =
                    g_restore_primary_metadata_dir + "/metadata_" +
                    std::to_string(task.file_idx) + ".bin";
                PrimaryMetadataReader metadata_reader(metadata_path);
                const std::string part_path =
                    parent_dir + "/restored_bam_parts/part_" +
                    std::to_string(task.file_idx) + ".bam";
                samFile* part_output = sam_open(part_path.c_str(), "wb");
                bool bam_part_ok = metadata_reader.valid() &&
                                   metadata_reader.block_count() == total_pairs &&
                                   part_output &&
                                   sam_hdr_write(part_output, restore_bam_header) >= 0;
                for (size_t i = 0; i < total_pairs && bam_part_ok; ++i) {
                    BamMetadataBlock metadata_block;
                    if (!metadata_reader.next(metadata_block) ||
                        metadata_block.is_paired != (blocks[i].is_paired != 0)) {
                        bam_part_ok = false;
                        break;
                    }
                    bam1_t* read1 = rebuild_bam_record(
                        metadata_block.read1, qnames[i],
                        reverse_pairs[i].first, quality_pairs[i].first);
                    if (!read1 ||
                        sam_write1(part_output, restore_bam_header, read1) < 0) {
                        bam_destroy1(read1);
                        bam_part_ok = false;
                        break;
                    }
                    bam_destroy1(read1);

                    if (metadata_block.is_paired) {
                        bam1_t* read2 = rebuild_bam_record(
                            metadata_block.read2, qnames[i],
                            reverse_pairs[i].second, quality_pairs[i].second);
                        if (!read2 ||
                            sam_write1(part_output, restore_bam_header, read2) < 0) {
                            bam_destroy1(read2);
                            bam_part_ok = false;
                            break;
                        }
                        bam_destroy1(read2);
                    }
                }
                if (part_output) sam_close(part_output);
                if (!bam_part_ok) {
                    err("BAM 分片重建失败: " + part_path, task.file_idx);
                    bam_restore_failed.store(true, std::memory_order_relaxed);
                }
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

    unsigned int requested_threads =
        share->o->thread_n > 0 ? (unsigned int)share->o->thread_n : 1u;

    unsigned int thread_count = std::max(1u, std::min(requested_threads, (unsigned int)total_tasks));

    std::cout << "[MAIN] requested_threads=" << requested_threads
              << ", thread_count=" << thread_count << std::endl;

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
    std::string output_prefix = share->o->fastq_output_dir
        ? std::string(share->o->fastq_output_dir)
        : parent_dir_final + "/merged_all";
    std::string merged_R1_path   = output_prefix + "_R1.fastq";
    std::string merged_R2_path   = output_prefix + "_R2.fastq";

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

    for (int id = 0; id <= file_num; ++id) {
        append(merged_r1, final_dir + "/output_" + std::to_string(id) + "_R1.fastq");
        if (g_input_is_paired) {
            append(merged_r2, final_dir + "/output_" + std::to_string(id) + "_R2.fastq");
        }
    }

    std::string error_fastq_path = parent_dir_final + "/error.fastq";
    uint64_t next_error_generated_id =
        (!partition_read_offsets.empty() ? partition_read_offsets.back() : 0) + 1;
    if (g_input_is_paired) {
        append_error_fastq_interleaved(error_fastq_path, merged_r1, merged_r2,
                                       share->o->preserve_qname == 0,
                                       next_error_generated_id);
    } else {
        append_error_fastq_single(error_fastq_path, merged_r1,
                                  share->o->preserve_qname == 0,
                                  next_error_generated_id);
    }

    merged_r1.close();
    if (g_input_is_paired) merged_r2.close();

    std::cout << "[MERGE] finish merge fastq" << std::endl;
    if (g_restore_bam_enabled &&
        !bam_restore_failed.load(std::memory_order_relaxed) &&
        !merge_restored_bam(file_num, parent_dir_final, restore_bam_header)) {
        std::cerr << "[BAM-RESTORE] failed to create restored BAM" << std::endl;
    }
    bam_hdr_destroy(restore_bam_header);
    std::cout << "[MAIN] go_to_decompress_pos finished, elapsed="
              << (cputime() - start_time) << " sec" << std::endl;
}

int DECOMPRESS_MAIN::run(int argc, char* argv[])
{
    std::cout << "work decompress" << std::endl;
    double cpu_time = cputime();
    double real_time1 = realtime();
    CLASSIFY_SHARE_DATA* share = NULL;
    share = new CLASSIFY_SHARE_DATA();
    g_share = share;
    share->o = new OL_PAR();
    if (share->o->get_decompress_option(argc, argv) != 0) {
        delete share->o;
        delete share;
        return 2;
    }
    const std::string archive_path = share->o->archive_path;
    const std::string parent_dir = archive_path + ".unpacked";
    const fs::path archive_parent = fs::absolute(fs::path(archive_path)).parent_path();
    const std::string archive_stem = fs::path(archive_path).stem().string();
    const std::string fastq_output_prefix =
        (archive_parent / archive_stem).string();
    share->o->fastq_output_dir = strdup(fastq_output_prefix.c_str());
    if (!extract_xzip_archive(archive_path, parent_dir)) {
        delete share->o; delete share;
        return 1;
    }
    if (!read_xzip_manifest(parent_dir, share->o)) {
        std::cerr << "[FATAL] invalid or missing xzip_manifest.txt" << std::endl;
        delete share->o; delete share;
        return 1;
    }
    g_global_o = share->o;
    g_restore_bam_enabled = false;
    if (share->o->restore_bam) {
        if (!share->o->bam_info_full || !share->o->preserve_qname) {
            std::cerr << "[FATAL] Cannot restore BAM: archive lacks full BAM metadata or preserved qnames"
                      << std::endl;
            delete share->o; delete share;
            return 1;
        }
        g_restore_bam_enabled = prepare_bam_restore_inputs(archive_path, parent_dir);
        if (!g_restore_bam_enabled) {
            delete share->o; delete share;
            return 1;
        }
    }
    if (g_restore_bam_enabled) {
        std::cout << "[BAM-RESTORE] enabled; output will be "
                  << g_restore_bam_output_path << std::endl;
    }
    share->o->pos_dir = strdup((parent_dir + "/final_store_with_name").c_str());
    share->o->webp_dir = strdup((parent_dir + "/store_2_zero_lossless").c_str());
    if (share->o->preserve_qname && !directory_exists(parent_dir + "/qnames")) {
        std::cerr << "[FATAL] archive manifest says qnames are preserved, but qnames/ is missing"
                  << std::endl;
        delete share->o; delete share;
        return 1;
    }
    g_input_is_paired = (share->o->paired_end != 0);
    char output_dir[512];
    if (create_output_dir_c(share->o->webp_dir, output_dir, sizeof(output_dir)) == 0) {
        printf("后续操作可用路径：%s\n", output_dir);
    }
    g_global_read_length = share->o->read_length;
    std::vector<uint32_t> pre_window_id_vector;  // 改名
    std::string pre_window_id_txt_path = parent_dir + "/pre_window_id.txt";
    std::ifstream pre_window_txt_file(pre_window_id_txt_path);
    if (!pre_window_txt_file.is_open()) {
        std::cerr << "错误：无法读取 " << pre_window_id_txt_path << std::endl;
        delete share->o; delete share;
        return 1;
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
        delete share->o; delete share;
        return 1;
    }
    std::cout << "成功读取 " << pre_window_id_vector.size() << " 个 pre_window_id。" << std::endl;
    genome = read_fasta(std::string(share->o->fasta_path));
    chr_bg_wb_ID = build_chr_bg_wb_ID();
    generate_window_info(chr_bg_wb_ID, share->o->thread_n);
    const QualityCodec quality_codec = read_quality_codec(parent_dir, false);
    if (quality_codec == QualityCodec::ZstdLossless) {
        std::cout << "[Quality] zstd lossless mode: use quality_score directory" << std::endl;
        if (!directory_exists(parent_dir + "/quality_score")) {
            std::cerr << "[FATAL] .xzip 中缺少 quality_score 目录" << std::endl;
            delete share->o; delete share;
            return 1;
        }
    } else {
        std::cout << "[Quality] " << quality_codec_name(quality_codec)
                  << " mode: reconstruct quality scores from WebP" << std::endl;
        int webp_rc = webp_reconstructor_40_main(
                share->o->webp_dir,
                output_dir,
                share->o->thread_n,
                share->o->read_length,
                chr_bg_wb_ID[genome.size()] / 10000
        );
        if (webp_rc != 0) {
            std::cerr << "[FATAL] WebP quality score reconstruction failed" << std::endl;
            delete share->o; delete share;
            return 1;
        }
    }
    share->o->webp_lossless = (quality_codec == QualityCodec::ZstdLossless) ? 1 : 0;
    go_to_decompress_pos(pre_window_id_vector, share);
    if (!share->o->keep_temp) {
        std::error_code cleanup_ec;
        fs::remove_all(parent_dir, cleanup_ec);
        if (cleanup_ec) {
            std::cerr << "[XZIP] failed to remove decompression temp directory: "
                      << cleanup_ec.message() << std::endl;
        }
    }
    std::cout << "run decompress successfully" << std::endl;
    fprintf(stderr, "Total CPU time: %.2f s\n", cputime() - cpu_time);
    fprintf(stderr, "Total REAL time: %.2f s\n", realtime() - real_time1);
    delete share->o;
    delete share;
    return 0;
}
