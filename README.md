# XZIP

XZIP is a reference-based compressor and decompressor for sequencing reads. The current implementation focuses on BAM input and reconstructs FASTQ output by storing each read as a compact combination of reference-window position, mismatch bitmap, mismatch bases, read name, orientation flag, and quality-score data.

The core implementation is in `src/BWT_aln.cpp`. The `Debug/` directory is an Eclipse/CDT-style debug build directory and contains generated makefiles, object files, logs, and test artifacts.

## Overview

XZIP compresses reads by using a reference FASTA as the reconstruction backbone:

1. Load the input BAM and reference FASTA.
2. Split each chromosome into fixed windows whose size is `read_length`.
3. For each read, calculate:
   - `window_id`: the reference window containing the read.
   - `hap_offset`: the effective offset inside that window, adjusted for leading soft/hard clipping.
   - `diff_seq`: a 0/1 bitmap where `0` means the read base equals the reference base and `1` means it differs.
   - `diff_base`: the actual bases at positions marked `1`.
4. Store read-pair position deltas, mismatch data, names, flags, and quality scores in separate files.
5. During decompression, rebuild each read from `window_id + hap_offset + diff_seq + diff_base + FASTA`, then restore quality strings, read names, and reverse-complement orientation.

This design avoids storing full read sequences when most bases match the reference.

## Main Commands

The executable dispatches subcommands in `src/main.cpp`:

```bash
./DNA_online_off_gz_yx bwt_aln [options] <idx_dir> <read_bam> <work_dir> <fasta_path> [rg]
./DNA_online_off_gz_yx decompress [options] <idx_dir> <pos_dir> <webp_dir> <fasta_path>
```

Depending on your build target, the executable name may differ. In the checked-in debug makefile, the linker target is currently `DNA_online_off_gz_yx`.

## Build

The debug build directory contains an auto-generated makefile:

```bash
cd Debug
make
```

The build links against zlib, pthread, WebP, and sharpyuv. The current debug makefile contains local WebP paths:

```make
-I/home/user/yexiang/libwebp/include
-L/home/user/yexiang/libwebp/lib
```

Adjust these paths if WebP is installed somewhere else.

## Compression

Example:

```bash
./DNA_online_off_gz_yx bwt_aln \
  -t 8 \
  -n 150 \
  -W 1 \
  -Z 19 \
  --paired_end 1 \
  /path/to/index_dir \
  /path/to/input.bam \
  /path/to/work_dir \
  /path/to/reference.fa
```

Important options:

| Option | Meaning | Default |
| --- | --- | --- |
| `-t`, `--thread` | Worker thread count | `4` |
| `-n`, `--read_length`, `--rlen` | Read length used as the reference-window size | `100` for compression |
| `-e`, `--paired_end` | `1` for paired-end BAM, `0` for single-end BAM | `1` |
| `-W`, `--webp_lossless`, `--quality_lossless` | `1` stores quality scores losslessly with zstd; `0` uses WebP lossy mode | `1` |
| `-Q`, `--webp_quality` | WebP quality when `-W 0` | `80` |
| `-M`, `--webp_method` | WebP method when `-W 0` | `6` |
| `-x`, `--webp_width` | WebP tile width, rounded to a multiple of read length | `4096` |
| `-y`, `--webp_height` | WebP tile height | `4096` |
| `-Z`, `--zstd_level` | zstd compression level | `19` |

The compression run creates subdirectories under `<work_dir>` automatically.

## Compression Output Layout

The main output files and directories are:

| Path under `work_dir` | Purpose |
| --- | --- |
| `tmp/test.<id>.bin` | Temporary partitioned `Compress_block` data grouped by `window_id / 10000` |
| `final_store_with_name/store.<id>.bin` | Huffman-coded position metadata, including window deltas and offsets |
| `diff_seq/diff_<id>.bin` | Packed 0/1 mismatch bitmaps |
| `diff_base/diff_<id>.bin` | Huffman-coded mismatch bases |
| `quality_score/quality_score.<id>.bin` | Raw quality blocks before optional zstd/WebP handling |
| `byte_flags/flags_<id>.bin` | Packed reverse-complement flag bits |
| `qnames/qnames_<id>.txt` | Read names |
| `pre_window_id.txt` | Previous-window anchors used to decode window deltas |
| `error.fastq` | Reads that could not be represented by the main path |
| `store_2_zero_lossless/` | WebP quality-score output when using lossy quality mode |

After compression, several directories/files may also be packaged as `.tar.zst`.

## Decompression

Example for lossless quality-score mode:

```bash
./DNA_online_off_gz_yx decompress \
  -t 8 \
  -n 150 \
  -W 1 \
  --paired_end 1 \
  /path/to/index_dir \
  /path/to/work_dir/final_store_with_name \
  /path/to/work_dir/store_2_zero_lossless \
  /path/to/reference.fa
```

The decompression options `read_length`, `paired_end`, and `webp_lossless` must match the compression run.

Decompression reads `pre_window_id.txt`, restores the reference windows from the FASTA, decodes each `store.<id>.bin`, applies `diff_seq` and `diff_base`, restores quality strings, read names, and reverse flags, and writes:

```text
work_dir/final_fastq/output_<id>_R1.fastq
work_dir/final_fastq/output_<id>_R2.fastq
work_dir/merged_all_R1.fastq
work_dir/merged_all_R2.fastq
```

For single-end mode, only the R1 output is produced. Records from `error.fastq` are appended to the merged FASTQ output.

## Data Model

The key compressed unit is `Compress_block`:

- R1 fields: `window_id`, `hap_offset`, `qname`, `flags_1`, `quality_score1`, `diff_seq1`, `diff_base1`.
- R2 fields: `window_id2`, `hap_offset2`, `flags_2`, `quality_score2`, `diff_seq2`, `diff_base2`.
- `is_paired` marks whether the block stores one read or a read pair.

Position metadata is converted to an A/C/G/T alphabet before Huffman coding:

- `window_dev`: current R1 window minus previous R1 window in the partition stream.
- `hap_offset`: R1 reference offset inside the window.
- `pair_window_dev`: absolute difference between R1 and R2 window IDs.
- `hap_offset2`: R2 reference offset.
- `add_flag`: whether R2 window is greater than or equal to R1 window.

The mismatch bitmap is stored separately because it has fixed length and packs efficiently as bits.

## Notes and Current Limitations

- Input is BAM-oriented. FASTQ reconstruction is the decompression output.
- The reference FASTA used for decompression must be the same as the one used for compression.
- `read_length` must match the actual read length expected by the run.
- In lossless quality mode, quality scores are packaged with zstd. In lossy mode, quality scores are represented through WebP-based image encoding.
- `idx_dir` remains part of the command-line interface, but the current reference-difference path primarily relies on BAM coordinates and the FASTA sequence.
- `Debug/` contains generated build outputs. Treat `src/` as the primary source tree for development.

