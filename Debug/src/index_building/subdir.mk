################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/index_building/haplotype_online.cpp 

CPP_DEPS += \
./src/index_building/haplotype_online.d 

OBJS += \
./src/index_building/haplotype_online.o 


# Each subdirectory must supply rules for building sources it contributes
src/index_building/%.o: ../src/index_building/%.cpp src/index_building/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I./htslib/htslib -I../ -I./htslib -I./htslib/cram -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-index_building

clean-src-2f-index_building:
	-$(RM) ./src/index_building/haplotype_online.d ./src/index_building/haplotype_online.o

.PHONY: clean-src-2f-index_building

