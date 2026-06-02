################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BWT_idx/bwt.cpp \
../src/BWT_idx/var_map.cpp 

CPP_DEPS += \
./src/BWT_idx/bwt.d \
./src/BWT_idx/var_map.d 

OBJS += \
./src/BWT_idx/bwt.o \
./src/BWT_idx/var_map.o 


# Each subdirectory must supply rules for building sources it contributes
src/BWT_idx/%.o: ../src/BWT_idx/%.cpp src/BWT_idx/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I./htslib/htslib -I../ -I./htslib -I./htslib/cram -O0 -Wall -g  -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-BWT_idx

clean-src-2f-BWT_idx:
	-$(RM) ./src/BWT_idx/bwt.d ./src/BWT_idx/bwt.o ./src/BWT_idx/var_map.d ./src/BWT_idx/var_map.o

.PHONY: clean-src-2f-BWT_idx

