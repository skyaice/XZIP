# XZIP

XZIP is a reference-based compressor for aligned sequencing reads. It encodes
read positions, sequence differences, pairing information, orientation flags,
and quality scores in a single `.xzip` archive. Decompression reconstructs
FASTQ files with the same reference genome used during compression.

## Requirements

- A C++17 compiler and a C compiler
- zlib
- libwebp and libsharpyuv
- `tar` and `zstd` on `PATH`

## Build

On macOS with Homebrew:

```bash
brew install webp zstd
make -j
```

For WebP installed outside the standard search path:

```bash
make -j WEBP_PREFIX=/path/to/libwebp
```

The executable is written to `build/xzip`.

## Compress

```bash
build/xzip compress [options] input.bam output_dir reference.fa
```

Example for paired-end reads with lossless quality scores:

```bash
build/xzip compress -t 8 -n 150 -W 1 \
  input.bam sample_xzip reference.fa
```

The final archive is `sample_xzip/sample_xzip.xzip`. Temporary compression
files are removed after the archive is created successfully.

Common options:

| Option | Description | Default |
| --- | --- | --- |
| `-t`, `--thread` | Worker threads | `4` |
| `-n`, `--read_length` | Read length | `100` |
| `-e`, `--paired_end` | `1` for paired-end, `0` for single-end | `1` |
| `-W`, `--webp_lossless` | `1` for lossless quality scores, `0` for lossy WebP | `1` |
| `-Q`, `--webp_quality` | WebP quality in lossy mode | `80` |
| `-M`, `--webp_method` | WebP compression method (`0`-`6`) | `6` |
| `-U`, `--preserve_qname` | Preserve original read names | off |
| `-I`, `--bam_info full` | Store metadata required for optional BAM restoration | `none` |

`bwt_aln` remains available as an alias for `compress` for compatibility with
older scripts.

## Decompress

```bash
build/xzip decompress [options] archive.xzip reference.fa
```

Example:

```bash
build/xzip decompress -t 8 sample_xzip/sample_xzip.xzip reference.fa
```

Read length, paired-end layout, quality codec, and read-name mode are read from
the archive manifest. Output files are written next to the archive:

```text
sample_xzip_R1.fastq
sample_xzip_R2.fastq
```

For single-end data, only the R1 file is produced. Archives created with
`--bam_info full` can also restore BAM output:

```bash
build/xzip decompress --restore_bam archive.xzip reference.fa
```

Use `--keep_temp` to retain the unpacked working directory for diagnostics.

## Archive contents

XZIP stores position deltas, mismatch bitmaps and bases, orientation flags,
quality-score data, and a manifest describing the decoding parameters. Read
names are optional. BAM restoration metadata is included only when requested.

The reference FASTA is not embedded in the archive and must match the reference
used for compression.

## Test

```bash
make check
```
