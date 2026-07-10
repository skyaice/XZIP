SHELL := /bin/sh

TARGET := build/xzip
BUILD_DIR := build/obj
WEBP_PREFIX ?= $(shell brew --prefix webp 2>/dev/null)

CXX ?= g++
CC ?= gcc

CPPFLAGS := -Isrc -Isrc/htslib -Isrc/htslib/htslib -Isrc/htslib/cram \
	-I/usr/local/include -I/opt/homebrew/include
CXXFLAGS := -std=c++17 -O3 -DNDEBUG -Wall -Wextra -MMD -MP
CFLAGS := -O3 -DNDEBUG -Wall -MMD -MP
LDFLAGS := -L/usr/local/lib -L/opt/homebrew/lib
LDLIBS := -lm -lpthread -lz -lwebp -lsharpyuv

ifneq ($(strip $(WEBP_PREFIX)),)
CPPFLAGS += -I$(WEBP_PREFIX)/include
LDFLAGS += -L$(WEBP_PREFIX)/lib
endif

CPP_SOURCES := \
	src/main.cpp \
	src/BWT_aln.cpp \
	src/webp_tool/compress_40_quality_score.cpp \
	src/webp_tool/webp_reconstructor_40.cpp

C_SOURCES := \
	src/clib/bam_file.c \
	src/clib/utils.c \
	$(wildcard src/htslib/*.c) \
	$(wildcard src/htslib/cram/*.c)

CPP_OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(CPP_SOURCES))
C_OBJECTS := $(patsubst %.c,$(BUILD_DIR)/%.o,$(C_SOURCES))
OBJECTS := $(CPP_OBJECTS) $(C_OBJECTS)
DEPS := $(OBJECTS:.o=.d)

.PHONY: all check clean

all: check-deps $(TARGET)

check-deps:
	@command -v zstd >/dev/null 2>&1 || { echo "Missing zstd executable"; exit 1; }
	@printf '#include <webp/decode.h>\n' | \
		$(CXX) $(CPPFLAGS) -x c++ -E - >/dev/null 2>&1 || \
		{ echo "Missing libwebp headers; set WEBP_PREFIX=/path/to/libwebp"; exit 1; }

$(TARGET): $(OBJECTS)
	@mkdir -p $(dir $@)
	$(CXX) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

check: all
	sh tests/smoke_cli.sh

clean:
	rm -rf build

-include $(DEPS)
