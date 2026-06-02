################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CPPLIB/Assembler/assembler.cpp 

CPP_DEPS += \
./src/CPPLIB/Assembler/assembler.d 

OBJS += \
./src/CPPLIB/Assembler/assembler.o 


# Each subdirectory must supply rules for building sources it contributes
src/CPPLIB/Assembler/%.o: ../src/CPPLIB/Assembler/%.cpp src/CPPLIB/Assembler/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I./htslib/htslib -I../ -I./htslib -I./htslib/cram -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-CPPLIB-2f-Assembler

clean-src-2f-CPPLIB-2f-Assembler:
	-$(RM) ./src/CPPLIB/Assembler/assembler.d ./src/CPPLIB/Assembler/assembler.o

.PHONY: clean-src-2f-CPPLIB-2f-Assembler

