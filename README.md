# XZIP

XZIP is a Linux-first, reference-based compressor for aligned sequencing reads.
It turns an aligned BAM plus its reference FASTA into one portable `.xzip`
archive and reconstructs FASTQ files from that archive. Instead of storing every
read sequence again, XZIP stores its reference position and only the bases that
differ from the reference, together with pairing, orientation, name (optional),
and quality information.

The reference FASTA is deliberately not embedded in the archive. Keep the exact
same reference available for decompression.

## Install on Linux

### 1. Install prerequisites

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

On RHEL, Rocky Linux, AlmaLinux, or CentOS, install the equivalent packages with
the package manager available on the server. XZIP requires a compiler with
C++17 support. Older CentOS installations often provide an older compiler by
default; activate a newer GCC toolset before building.

If you do not have administrator privileges, ask the server administrator for
the following development/runtime dependencies: a C++17 compiler, zlib,
libwebp, zstd, and tar. A `zstd` executable alone is not enough: the libwebp
headers (for example `webp/decode.h`) must also be installed.

On macOS with Homebrew:

```bash
brew install webp zstd
```

### 2. Clone, check, and build

```bash
git clone https://github.com/skyaice/XZIP.git
cd XZIP
make doctor
make -j"$(nproc)"
```

The executable is written to `build/xzip`. The Makefile uses `pkg-config` to
find libwebp on Linux. If libwebp is installed in a custom location, provide
its prefix explicitly:

```bash
make -j WEBP_PREFIX=/path/to/libwebp
```

Verify the build before using real data:

```bash
make check
build/xzip --help
```

To make `xzip` available from any directory without administrator access:

```bash
make install PREFIX="${HOME}/.local"
export PATH="${HOME}/.local/bin:${PATH}"
xzip --help
```

Add the `export PATH=...` line to your shell profile if it is not already
present. System administrators can use the default `/usr/local` prefix with
appropriate permissions.

The bundled `src/htslib` directory includes the complete upstream configuration
and packaging files needed to rebuild or inspect the vendored dependency on a
Linux system. XZIP itself uses the checked-in `src/htslib/config.h` for its
normal build.

## Five-minute quick start

Compress a paired-end BAM with 150 bp reads using eight threads:

```bash
xzip compress -t 8 -n 150 -W 1 \
  input.bam sample_xzip reference.fa
```

This creates `sample_xzip/sample_xzip.xzip`. After successful packaging,
intermediate files are removed from `sample_xzip`.

Reconstruct the FASTQ files:

```bash
xzip decompress -t 8 \
  sample_xzip/sample_xzip.xzip reference.fa
```

The output is written beside the archive:

```text
sample_xzip_R1.fastq
sample_xzip_R2.fastq
```

For single-end input, only `sample_xzip_R1.fastq` is produced.

If you skipped `make install`, replace `xzip` in every example with
`build/xzip`.

## How the compression works

For each primary BAM read, XZIP:

1. Locates the read in a fixed-size reference window determined by the read
   length.
2. Stores the window delta and the read's offset inside that window instead of
   an absolute sequence copy.
3. Compares the read with the reference, bit-packs a match/mismatch bitmap, and
   Huffman-encodes the bases at mismatch positions.
4. Stores pairing distance, reverse-complement flags, optional read names, and
   quality scores in separate streams so each stream can use an appropriate
   codec.
5. Packages the manifest and streams as a tar archive compressed with zstd.

Decompression reverses those steps: it restores the reference window, applies
the mismatch bitmap and bases, restores orientation and quality data, and writes
R1/R2 FASTQ files. Reads that cannot use the main reference-difference path are
preserved through a fallback stream rather than silently discarded.

## Compression command

```text
xzip compress [options] <input.bam> <output_dir> <reference.fa> [read_group]
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

## Decompression command

```text
xzip decompress [options] <archive.xzip> <reference.fa>
```

The archive manifest supplies the read length, pairing layout, quality codec,
and read-name mode. Use `--keep_temp` to retain the unpacked working directory
for inspection.

### Optional BAM restoration

To request a restored BAM, create the archive with full metadata:

```bash
xzip compress -I full -n 150 input.bam sample_xzip reference.fa
```

Then run:

```bash
xzip decompress --restore_bam sample_xzip/sample_xzip.xzip reference.fa
```

Full metadata preserves the information required by the restoration workflow
and enables read-name preservation automatically.

## Operational checklist

Before a production run, confirm all of the following:

- The BAM was aligned against the supplied FASTA.
- `-n` exactly matches the read length in the BAM.
- `-e 1` is used for paired-end data and `-e 0` for single-end data.
- The same, unchanged FASTA will be retained for decompression.
- There is enough temporary disk space for the working directory and final
  archive.
- Lossy quality mode (`-W 0`) is used only when approximate quality scores are
  acceptable.

## Validation and troubleshooting

Build and run the command-line smoke test with:

```bash
make check
```

For a new dataset, first test a small BAM and compare the reconstructed FASTQ
record count, sequences, and quality values with the expected result.

Common build failures:

| Message | Resolution |
| --- | --- |
| `Missing zstd executable` | Install the zstd command-line package and ensure it is on `PATH`. |
| `Missing libwebp headers` | Install the libwebp development package, or pass `WEBP_PREFIX=/path/to/libwebp`. |
| C++17 compilation errors | Activate a newer GCC/Clang toolchain, then run `make clean && make`. |
| Permission denied during `make install` | Use `PREFIX="${HOME}/.local"` or ask an administrator to install system-wide. |

Run `make doctor` whenever the project is moved to a new server. It reports the
compiler and dependency paths used by the build.
