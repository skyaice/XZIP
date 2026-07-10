#include "webp_reconstructor_40.h"

#include <webp/decode.h>

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;

namespace {

int block_number(const fs::path& path) {
    const std::string name = path.stem().string();
    const size_t pos = name.find_last_of('_');
    if (pos == std::string::npos) return -1;
    try {
        return std::stoi(name.substr(pos + 1));
    } catch (...) {
        return -1;
    }
}

bool reconstruct_one_id(int id,
                        const std::string& webp_base_dir,
                        const std::string& output_dir,
                        int read_length,
                        std::string& error) {
    const fs::path input_dir =
        fs::path(webp_base_dir) / ("combined_output." + std::to_string(id));
    if (!fs::is_directory(input_dir)) return true;

    size_t sequence_count = 0;
    {
        std::ifstream count_file(input_dir / "sequence_count.txt");
        if (!(count_file >> sequence_count)) {
            error = "missing sequence_count.txt in " + input_dir.string();
            return false;
        }
    }

    std::vector<fs::path> webp_files;
    for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".webp") {
            webp_files.push_back(entry.path());
        }
    }
    std::sort(webp_files.begin(), webp_files.end(),
              [](const fs::path& a, const fs::path& b) {
                  return block_number(a) < block_number(b);
              });

    std::vector<std::string> qualities;
    qualities.reserve(sequence_count);
    for (const fs::path& path : webp_files) {
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        if (!in) {
            error = "cannot open " + path.string();
            return false;
        }
        const std::streamsize size = in.tellg();
        in.seekg(0);
        std::vector<uint8_t> encoded(static_cast<size_t>(size));
        if (size <= 0 || !in.read(reinterpret_cast<char*>(encoded.data()), size)) {
            error = "cannot read " + path.string();
            return false;
        }

        int width = 0;
        int height = 0;
        uint8_t* rgba = WebPDecodeRGBA(encoded.data(), encoded.size(), &width, &height);
        if (!rgba || width <= 0 || height <= 0 || width % read_length != 0) {
            if (rgba) WebPFree(rgba);
            error = "invalid WebP dimensions in " + path.string();
            return false;
        }

        const int sequences_per_row = width / read_length;
        for (int y = 0; y < height && qualities.size() < sequence_count; ++y) {
            for (int seq = 0; seq < sequences_per_row && qualities.size() < sequence_count; ++seq) {
                std::string quality(static_cast<size_t>(read_length), '\0');
                const int x_start = seq * read_length;
                for (int x = 0; x < read_length; ++x) {
                    quality[static_cast<size_t>(x)] =
                        static_cast<char>(rgba[(static_cast<size_t>(y) * width + x_start + x) * 4]);
                }
                qualities.push_back(std::move(quality));
            }
        }
        WebPFree(rgba);
    }

    if (qualities.size() != sequence_count) {
        error = "sequence count mismatch for ID " + std::to_string(id);
        return false;
    }

    fs::create_directories(output_dir);
    const fs::path output_path =
        fs::path(output_dir) / ("quality_score." + std::to_string(id) + ".bin");
    std::ofstream out(output_path, std::ios::binary | std::ios::trunc);
    if (!out) {
        error = "cannot create " + output_path.string();
        return false;
    }
    out.write(reinterpret_cast<const char*>(&sequence_count), sizeof(sequence_count));
    for (const std::string& quality : qualities) {
        const size_t len = quality.size();
        out.write(reinterpret_cast<const char*>(&len), sizeof(len));
        out.write(quality.data(), static_cast<std::streamsize>(quality.size()));
    }
    if (!out.good()) {
        error = "cannot write " + output_path.string();
        return false;
    }
    return true;
}

}  // namespace

int webp_reconstructor_40_main(char* webp_dir,
                               char* output_dir,
                               int thread_n,
                               int read_length,
                               int total) {
    if (!webp_dir || !output_dir || read_length <= 0 || total < 0) return 1;

    const int workers = std::max(1, thread_n);
    std::atomic<int> next_id{0};
    std::atomic<bool> failed{false};
    std::mutex error_mutex;
    std::string first_error;
    std::vector<std::thread> threads;

    for (int i = 0; i < workers; ++i) {
        threads.emplace_back([&] {
            while (!failed.load()) {
                const int id = next_id.fetch_add(1);
                if (id > total) break;
                std::string error;
                if (!reconstruct_one_id(id, webp_dir, output_dir, read_length, error)) {
                    failed = true;
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (first_error.empty()) first_error = std::move(error);
                }
            }
        });
    }
    for (auto& thread : threads) thread.join();

    if (failed) {
        std::cerr << "[WebP reconstruction] " << first_error << std::endl;
        return 1;
    }
    return 0;
}
