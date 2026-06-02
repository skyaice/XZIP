################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/clib/kswlib/kalloc.c \
../src/clib/kswlib/ksw2_dispatch.c \
../src/clib/kswlib/ksw2_extd2_sse.c 

C_DEPS += \
./src/clib/kswlib/kalloc.d \
./src/clib/kswlib/ksw2_dispatch.d \
./src/clib/kswlib/ksw2_extd2_sse.d 

OBJS += \
./src/clib/kswlib/kalloc.o \
./src/clib/kswlib/ksw2_dispatch.o \
./src/clib/kswlib/ksw2_extd2_sse.o 


# Each subdirectory must supply rules for building sources it contributes
src/clib/kswlib/%.o: ../src/clib/kswlib/%.c src/clib/kswlib/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -lc -static -static-libstdc++  -static-libgcc -I"../src/htslib/" -I"../src" -I"../src/htslib/cram" -I"../src/htslib/htslib" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-clib-2f-kswlib

clean-src-2f-clib-2f-kswlib:
	-$(RM) ./src/clib/kswlib/kalloc.d ./src/clib/kswlib/kalloc.o ./src/clib/kswlib/ksw2_dispatch.d ./src/clib/kswlib/ksw2_dispatch.o ./src/clib/kswlib/ksw2_extd2_sse.d ./src/clib/kswlib/ksw2_extd2_sse.o

.PHONY: clean-src-2f-clib-2f-kswlib

