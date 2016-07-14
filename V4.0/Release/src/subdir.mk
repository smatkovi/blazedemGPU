################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/main.cpp 

OBJS += \
./src/main.o 

CPP_DEPS += \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-7.5/bin/nvcc -lineinfo -O3 -ftz true -v -Xcompiler -fPIC -Xptxas -dlcm=ca -gencode arch=compute_30,code=sm_30 -gencode arch=compute_52,code=sm_52 -m64 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-7.5/bin/nvcc -lineinfo -O3 -ftz true -v -Xcompiler -fPIC -Xptxas -dlcm=ca --compile -m64  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


