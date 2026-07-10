#include "BWT_aln.hpp"
#include "CPPLIB/tools.hpp"

#include <cstring>
#include <iostream>

namespace {

void print_usage(const char* program) {
    std::cerr
        << "XZIP reference-based read compressor\n\n"
        << "Usage:\n"
        << "  " << program << " compress [options] <input.bam> <work_dir> <reference.fa> [read_group]\n"
        << "  " << program << " decompress [options] <archive.xzip> <reference.fa>\n\n"
        << "Run either command with --help for its full option list.\n";
}

int run_compress(int argc, char** argv) {
    BWT_aln::BWT_CLASSIFY_MAIN compressor;
    compressor.init_run(argc, argv);
    return 0;
}

int run_decompress(int argc, char** argv) {
    BWT_aln::DECOMPRESS_MAIN decompressor;
    decompressor.run(argc, argv);
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    union EndianCheck { int value; char byte; } endian;
    endian.value = 1;
    xassert(endian.byte == 1, "XZIP requires a little-endian system");

    if (argc < 2 || std::strcmp(argv[1], "--help") == 0 ||
        std::strcmp(argv[1], "-h") == 0) {
        print_usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }
    if (std::strcmp(argv[1], "compress") == 0 ||
        std::strcmp(argv[1], "bwt_aln") == 0) {
        return run_compress(argc - 1, argv + 1);
    }
    if (std::strcmp(argv[1], "decompress") == 0) {
        return run_decompress(argc - 1, argv + 1);
    }

    std::cerr << "Unknown command: " << argv[1] << "\n\n";
    print_usage(argv[0]);
    return 1;
}
