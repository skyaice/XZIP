# XZIP

XZIP is a reference-based compressor for aligned sequencing reads. It packages
the information needed to reconstruct FASTQ files from a BAM file, including
read placement, sequence differences, pairing, orientation, and quality scores.
The reference FASTA is supplied again at decompression time and is not embedded
in the archive.

## Requirements

- C++17 and C compilers
- zlib
- libwebp development headers and library
- `tar` and `zstd` available on `PATH`

On Debian or Ubuntu:

```bash
sudo apt update
sudo apt install -y build-essential pkg-config zlib1g-dev libwebp-dev zstd tar
```

On Fedora:

```bash
sudo dnf install -y gcc-c++ gcc make pkgconf-pkg-config zlib-devel libwebp-devel zstd tar
```

On macOS with Homebrew:

```bash
brew install webp zstd
```

## Build

```bash
git clone https://github.com/skyaice/XZIP.git
cd XZIP
make -j
```

The executable is written to `build/xzip`. The Makefile uses `pkg-config` to
find libwebp on Linux. If libwebp is installed in a custom location, provide
its prefix explicitly:

```bash
make -j WEBP_PREFIX=/path/to/libwebp
```

The bundled `src/htslib` directory includes the complete upstream configuration
and packaging files needed to rebuild or inspect the vendored dependency on a
Linux system. XZIP itself uses the checked-in `src/htslib/config.h` for its
normal build.

## Quick start

Compress a paired-end BAM with 150 bp reads using eight threads:

```bash
build/xzip compress -t 8 -n 150 -W 1 \
  input.bam sample_xzip reference.fa
```

This creates `sample_xzip/sample_xzip.xzip`. After successful packaging,
intermediate files are removed from `sample_xzip`.

Reconstruct the FASTQ files:

```bash
build/xzip decompress -t 8 \
  sample_xzip/sample_xzip.xzip reference.fa
```

The output is written beside the archive:

```text
sample_xzip_R1.fastq
sample_xzip_R2.fastq
```

For single-end input, only `sample_xzip_R1.fastq` is produced.

## Compression

```text
build/xzip compress [options] <input.bam> <output_dir> <reference.fa> [read_group]
```

`input.bam` must be aligned against `reference.fa`. The read length given by
`-n` must match the input data. Keep the reference FASTA unchanged: the same
file is required to decompress the archive.

Common options:

| Option | Meaning | Default |
| --- | --- | --- |
| `-t`, `--thread` | Number of worker threads | `4` |
| `-n`, `--read_length` | Read length in bases | `100` |
| `-e`, `--paired_end` | `1` for paired-end, `0` for single-end input | `1` |
| `-W`, `--webp_lossless` | `1` for lossless quality handling, `0` for lossy WebP quality compression | `1` |
| `-Q`, `--webp_quality` | WebP quality when `-W 0` is selected | `80` |
| `-M`, `--webp_method` | WebP compression method (`0`-`6`) | `6` |
| `-U`, `--preserve_qname` | Store original read names | off |
| `-I`, `--bam_info full` | Store metadata for optional BAM restoration | `none` |

With the default lossless setting, XZIP uses lossless WebP for compatible
four-value quality data and lossless zstd packaging otherwise. Setting `-W 0`
enables lossy WebP quality compression; use it only when approximate quality
scores are acceptable.

`bwt_aln` remains an alias for `compress` for compatibility with earlier
scripts.

## Decompression

```text
build/xzip decompress [options] <archive.xzip> <reference.fa>
```

The archive manifest supplies the read length, pairing layout, quality codec,
and read-name mode. Use `--keep_temp` to retain the unpacked working directory
for inspection.

### Optional BAM restoration

To request a restored BAM, create the archive with full metadata:

```bash
build/xzip compress -I full -n 150 input.bam sample_xzip reference.fa
```

Then run:

```bash
build/xzip decompress --restore_bam sample_xzip/sample_xzip.xzip reference.fa
```

Full metadata preserves the information required by the restoration workflow
and enables read-name preservation automatically.

## Validation

Build and run the command-line smoke test with:

```bash
make check
```

For a new dataset, first test a small BAM and compare the reconstructed FASTQ
record count, sequences, and quality values with the expected result.
