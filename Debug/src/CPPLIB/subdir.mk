################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CPPLIB/tools.cpp 

CPP_DEPS += \
./src/CPPLIB/tools.d 

OBJS += \
./src/CPPLIB/tools.o 


# Each subdirectory must supply rules for building sources it contributes
src/CPPLIB/%.o: ../src/CPPLIB/%.cpp src/CPPLIB/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I./htslib/htslib -I../ -I./htslib -I./htslib/cram -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-CPPLIB

clean-src-2f-CPPLIB:
	-$(RM) ./src/CPPLIB/tools.d ./src/CPPLIB/tools.o

.PHONY: clean-src-2f-CPPLIB

