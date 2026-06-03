TARGET ?= DNA_online_off_gz_yx
BUILD_DIR ?= build

CXX ?= g++
CC ?= gcc
RM ?= rm -rf

CXXFLAGS ?= -O3 -g -Wall -std=c++17
CFLAGS ?= -O3 -g -Wall

ZSTD_CFLAGS ?= $(shell pkg-config --cflags libzstd 2>/dev/null)
ZSTD_LIBS ?= $(shell pkg-config --libs libzstd 2>/dev/null)
ifeq ($(strip $(ZSTD_LIBS)),)
ZSTD_LIBS := -lzstd
endif

WEBP_CFLAGS ?= $(shell pkg-config --cflags libwebp 2>/dev/null)
WEBP_LIBS ?= $(shell pkg-config --libs libwebp libsharpyuv 2>/dev/null)
ifeq ($(strip $(WEBP_LIBS)),)
WEBP_LIBS := -lwebp -lsharpyuv
endif

CPPFLAGS += \
	-Isrc \
	-Isrc/htslib \
	-Isrc/htslib/cram \
	-Isrc/htslib/htslib \
	$(ZSTD_CFLAGS) \
	$(WEBP_CFLAGS)

LDLIBS += -lm -lpthread -lz $(ZSTD_LIBS) $(WEBP_LIBS)

CPP_SRCS := \
	src/ReadHandler.cpp \
	src/BWT_aln.cpp \
	src/aln_RST_DEF.cpp \
	src/main.cpp \
	src/occ_RST_DEF.cpp \
	src/process3_var_call.cpp \
	src/BWT_idx/bwt.cpp \
	src/BWT_idx/var_map.cpp \
	src/CPPLIB/Assembler/assembler.cpp \
	src/CPPLIB/JsonObject/CJsonObject.cpp \
	src/CPPLIB/JsonObject/demo.cpp \
	src/CPPLIB/tools.cpp \
	src/index_building/haplotype_online.cpp \
	src/minimizer/minimizer.cpp \
	src/var_calling/GT_TYPE.cpp \
	src/var_calling/known_var_candidate_generator.cpp \
	src/var_calling/known_var_hap_counter.cpp \
	src/var_calling/novel_var_candidate_generator.cpp \
	src/var_calling/var_caller_all_type.cpp \
	src/webp_tool/compress_40_quality_score.cpp \
	src/webp_tool/webp_converter.cpp \
	src/webp_tool/webp_reconstructor_4.cpp \
	src/webp_tool/webp_reconstructor_40.cpp

C_SRCS := \
	src/CPPLIB/JsonObject/cJSON.c \
	src/bwa_idx/QSufSort.c \
	src/bwa_idx/bntseq.c \
	src/bwa_idx/bwt.c \
	src/bwa_idx/bwt_gen.c \
	src/bwa_idx/bwtindex.c \
	src/bwa_idx/is.c \
	src/bwa_idx/rle.c \
	src/bwa_idx/rope.c \
	src/clib/bam_file.c \
	src/clib/binarys_qsort.c \
	src/clib/kthread.c \
	src/clib/utils.c \
	src/clib/vcf_file.c \
	src/clib/kswlib/kalloc.c \
	src/clib/kswlib/ksw2_dispatch.c \
	src/clib/kswlib/ksw2_extd2_sse.c \
	src/htslib/bcf_sr_sort.c \
	src/htslib/bgzf.c \
	src/htslib/bgzip.c \
	src/htslib/errmod.c \
	src/htslib/faidx.c \
	src/htslib/hfile.c \
	src/htslib/hfile_gcs.c \
	src/htslib/hfile_net.c \
	src/htslib/hts.c \
	src/htslib/hts_os.c \
	src/htslib/htsfile.c \
	src/htslib/kfunc.c \
	src/htslib/knetfile.c \
	src/htslib/kstring.c \
	src/htslib/md5.c \
	src/htslib/multipart.c \
	src/htslib/probaln.c \
	src/htslib/realn.c \
	src/htslib/regidx.c \
	src/htslib/sam.c \
	src/htslib/synced_bcf_reader.c \
	src/htslib/tabix.c \
	src/htslib/tbx.c \
	src/htslib/textutils.c \
	src/htslib/thread_pool.c \
	src/htslib/vcf.c \
	src/htslib/vcf_sweep.c \
	src/htslib/vcfutils.c \
	src/htslib/cram/cram_codecs.c \
	src/htslib/cram/cram_decode.c \
	src/htslib/cram/cram_encode.c \
	src/htslib/cram/cram_external.c \
	src/htslib/cram/cram_index.c \
	src/htslib/cram/cram_io.c \
	src/htslib/cram/cram_samtools.c \
	src/htslib/cram/cram_stats.c \
	src/htslib/cram/files.c \
	src/htslib/cram/mFILE.c \
	src/htslib/cram/open_trace_file.c \
	src/htslib/cram/pooled_alloc.c \
	src/htslib/cram/rANS_static.c \
	src/htslib/cram/sam_header.c \
	src/htslib/cram/string_alloc.c

CPP_OBJS := $(patsubst src/%.cpp,$(BUILD_DIR)/src/%.o,$(CPP_SRCS))
C_OBJS := $(patsubst src/%.c,$(BUILD_DIR)/src/%.o,$(C_SRCS))
OBJS := $(CPP_OBJS) $(C_OBJS)
DEPS := $(OBJS:.o=.d)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDLIBS)

$(BUILD_DIR)/src/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c -o $@ $<

$(BUILD_DIR)/src/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c -o $@ $<

clean:
	$(RM) $(BUILD_DIR) $(TARGET)

-include $(DEPS)
