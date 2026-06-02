################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

CXXFLAGS += -std=c++17

ZSTD_CFLAGS ?= $(shell pkg-config --cflags libzstd 2>/dev/null)
ZSTD_LIBS ?= $(shell pkg-config --libs libzstd 2>/dev/null)
ifeq ($(strip $(ZSTD_LIBS)),)
ZSTD_LIBS := -lzstd
endif

WEBP_CONVERTER_SRC := ../webp_converter.cpp
COMPRESS_40_QUALITY_SCORE_SRC := ../compress_40_quality_score.cpp

# All of the sources participating in the build are defined here
-include sources.mk
-include src/var_calling/subdir.mk
-include src/minimizer/subdir.mk
-include src/index_building/subdir.mk
-include src/htslib/cram/subdir.mk
-include src/htslib/subdir.mk
-include src/clib/kswlib/subdir.mk
-include src/clib/subdir.mk
-include src/bwa_idx/subdir.mk
-include src/CPPLIB/JsonObject/subdir.mk
-include src/CPPLIB/Assembler/subdir.mk
-include src/CPPLIB/subdir.mk
-include src/BWT_idx/subdir.mk
-include src/subdir.mk
-include src/webp_tool/subdir.mk
-include src/genome_read/subdir.mk
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

OPTIONAL_TOOL_DEPS := \
$(wildcard ../makefile.defs) \
$(wildcard ../makefile.init) \
$(wildcard ../makefile.targets) \


BUILD_ARTIFACT_NAME := DNA_online_off_yx_1
BUILD_ARTIFACT_EXTENSION :=
BUILD_ARTIFACT_PREFIX :=
BUILD_ARTIFACT := $(BUILD_ARTIFACT_PREFIX)$(BUILD_ARTIFACT_NAME)$(if $(BUILD_ARTIFACT_EXTENSION),.$(BUILD_ARTIFACT_EXTENSION),)

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: main-build

# Main-build Target
main-build: DNA_online_off

# Tool invocations
DNA_online_off: $(OBJS) $(USER_OBJS) makefile $(OPTIONAL_TOOL_DEPS) 
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++ -o "DNA_online_off_gz_yx" $(OBJS) $(USER_OBJS) $(LIBS) -I/home/user/yexiang/libwebp/include -L/home/user/yexiang/libwebp/lib -lm -lpthread -lz $(ZSTD_LIBS) -lwebp -lsharpyuv
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) DNA_online_off
	-@echo ' '

.PHONY: all clean dependents main-build

-include ../makefile.targets
