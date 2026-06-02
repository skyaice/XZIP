################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CPPLIB/JsonObject/CJsonObject.cpp \
../src/CPPLIB/JsonObject/demo.cpp 

C_SRCS += \
../src/CPPLIB/JsonObject/cJSON.c 

CPP_DEPS += \
./src/CPPLIB/JsonObject/CJsonObject.d \
./src/CPPLIB/JsonObject/demo.d 

C_DEPS += \
./src/CPPLIB/JsonObject/cJSON.d 

OBJS += \
./src/CPPLIB/JsonObject/CJsonObject.o \
./src/CPPLIB/JsonObject/cJSON.o \
./src/CPPLIB/JsonObject/demo.o 


# Each subdirectory must supply rules for building sources it contributes
src/CPPLIB/JsonObject/%.o: ../src/CPPLIB/JsonObject/%.cpp src/CPPLIB/JsonObject/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I./htslib/htslib -I../ -I./htslib -I./htslib/cram -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/CPPLIB/JsonObject/%.o: ../src/CPPLIB/JsonObject/%.c src/CPPLIB/JsonObject/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -lc -static -static-libstdc++  -static-libgcc -I"../src/htslib/" -I"../src" -I"../src/htslib/cram" -I"../src/htslib/htslib" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src-2f-CPPLIB-2f-JsonObject

clean-src-2f-CPPLIB-2f-JsonObject:
	-$(RM) ./src/CPPLIB/JsonObject/CJsonObject.d ./src/CPPLIB/JsonObject/CJsonObject.o ./src/CPPLIB/JsonObject/cJSON.d ./src/CPPLIB/JsonObject/cJSON.o ./src/CPPLIB/JsonObject/demo.d ./src/CPPLIB/JsonObject/demo.o

.PHONY: clean-src-2f-CPPLIB-2f-JsonObject

