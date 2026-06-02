################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/bwa_idx/QSufSort.c \
../src/bwa_idx/bntseq.c \
../src/bwa_idx/bwt.c \
../src/bwa_idx/bwt_gen.c \
../src/bwa_idx/bwtindex.c \
../src/bwa_idx/is.c \
../src/bwa_idx/rle.c \
../src/bwa_idx/rope.c 

C_DEPS += \
./src/bwa_idx/QSufSort.d \
./src/bwa_idx/bntseq.d \
./src/bwa_idx/bwt.d \
./src/bwa_idx/bwt_gen.d \
./src/bwa_idx/bwtindex.d \
./src/bwa_idx/is.d \
./src/bwa_idx/rle.d \
./src/bwa_idx/rope.d 

OBJS += \
./src/bwa_idx/QSufSort.o \
./src/bwa_idx/bntseq.o \
./src/bwa_idx/bwt.o \
./src/bwa_idx/bwt_gen.o \
./src/bwa_idx/bwtindex.o \
./src/bwa_idx/is.o \
./src/bwa_idx/rle.o \
./src/bwa_idx/rope.o 


# Each subdirectory must supply rules for building sources it contributes
src/bwa_idx/%.o: ../src/bwa_idx/%.c src/bwa_idx/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I"../src/htslib/" -I"../src" -I"../src/htslib/cram" -I"../src/htslib/htslib" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-bwa_idx

clean-src-2f-bwa_idx:
	-$(RM) ./src/bwa_idx/QSufSort.d ./src/bwa_idx/QSufSort.o ./src/bwa_idx/bntseq.d ./src/bwa_idx/bntseq.o ./src/bwa_idx/bwt.d ./src/bwa_idx/bwt.o ./src/bwa_idx/bwt_gen.d ./src/bwa_idx/bwt_gen.o ./src/bwa_idx/bwtindex.d ./src/bwa_idx/bwtindex.o ./src/bwa_idx/is.d ./src/bwa_idx/is.o ./src/bwa_idx/rle.d ./src/bwa_idx/rle.o ./src/bwa_idx/rope.d ./src/bwa_idx/rope.o

.PHONY: clean-src-2f-bwa_idx

