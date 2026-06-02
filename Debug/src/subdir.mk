################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/ReadHandler.cpp \
../src/BWT_aln.cpp \
../src/aln_RST_DEF.cpp \
../src/main.cpp \
../src/occ_RST_DEF.cpp \
../src/process3_var_call.cpp 

CPP_DEPS += \
./src/ReadHandler.d\
./src/BWT_aln.d \
./src/aln_RST_DEF.d \
./src/main.d \
./src/occ_RST_DEF.d \
./src/process3_var_call.d 

OBJS += \
./src/ReadHandler.o\
./src/BWT_aln.o \
./src/aln_RST_DEF.o \
./src/main.o \
./src/occ_RST_DEF.o \
./src/process3_var_call.o 


# Each subdirectory must supply rules for building sources it contributes
# 关键：基于 Debug 目录（当前工作目录）定位所有头文件路径
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	# 1. -I../src/htslib：定位 htslib/vcf.h（实际路径：../src/htslib/htslib/vcf.h）
	# 2. -I../src/htslib/cram：定位 cram 相关头文件（若需）
	# 3. -I/home/user/yexiang/libwebp/include：定位 webp/decode.h
	# 4. -I../：定位项目自身的头文件（如 ReadHandler.hpp）
	g++ -I../src/htslib \
        -I../src/htslib/cram \
        -I/home/user/yexiang/libwebp/include \
        -I../ \
        -O0 -Wall -g -c -fmessage-length=0 -MMD -MP \
        -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/BWT_aln.d ./src/BWT_aln.o ./src/aln_RST_DEF.d ./src/aln_RST_DEF.o ./src/main.d ./src/main.o ./src/occ_RST_DEF.d ./src/occ_RST_DEF.o ./src/process3_var_call.d ./src/process3_var_call.o

.PHONY: clean-src