################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/minimizer/minimizer.cpp 

CPP_DEPS += \
./src/minimizer/minimizer.d 

OBJS += \
./src/minimizer/minimizer.o 


# Each subdirectory must supply rules for building sources it contributes
src/minimizer/%.o: ../src/minimizer/%.cpp src/minimizer/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I./htslib/htslib -I../ -I./htslib -I./htslib/cram -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-minimizer

clean-src-2f-minimizer:
	-$(RM) ./src/minimizer/minimizer.d ./src/minimizer/minimizer.o

.PHONY: clean-src-2f-minimizer

