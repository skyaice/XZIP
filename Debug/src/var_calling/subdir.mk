################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/var_calling/GT_TYPE.cpp \
../src/var_calling/known_var_candidate_generator.cpp \
../src/var_calling/known_var_hap_counter.cpp \
../src/var_calling/novel_var_candidate_generator.cpp \
../src/var_calling/var_caller_all_type.cpp 

CPP_DEPS += \
./src/var_calling/GT_TYPE.d \
./src/var_calling/known_var_candidate_generator.d \
./src/var_calling/known_var_hap_counter.d \
./src/var_calling/novel_var_candidate_generator.d \
./src/var_calling/var_caller_all_type.d 

OBJS += \
./src/var_calling/GT_TYPE.o \
./src/var_calling/known_var_candidate_generator.o \
./src/var_calling/known_var_hap_counter.o \
./src/var_calling/novel_var_candidate_generator.o \
./src/var_calling/var_caller_all_type.o 


# Each subdirectory must supply rules for building sources it contributes
src/var_calling/%.o: ../src/var_calling/%.cpp src/var_calling/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I./htslib/htslib -I../ -I./htslib -I./htslib/cram -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-var_calling

clean-src-2f-var_calling:
	-$(RM) ./src/var_calling/GT_TYPE.d ./src/var_calling/GT_TYPE.o ./src/var_calling/known_var_candidate_generator.d ./src/var_calling/known_var_candidate_generator.o ./src/var_calling/known_var_hap_counter.d ./src/var_calling/known_var_hap_counter.o ./src/var_calling/novel_var_candidate_generator.d ./src/var_calling/novel_var_candidate_generator.o ./src/var_calling/var_caller_all_type.d ./src/var_calling/var_caller_all_type.o

.PHONY: clean-src-2f-var_calling

