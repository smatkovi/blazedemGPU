################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/Device/DeviceInterface.cu 

CU_DEPS += \
./src/Device/DeviceInterface.d 

OBJS += \
./src/Device/DeviceInterface.o 


# Each subdirectory must supply rules for building sources it contributes
src/Device/%.o: ../src/Device/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-7.5/bin/nvcc -lineinfo -O3 -ftz true -v -Xcompiler -fPIC -Xptxas -dlcm=ca -gencode arch=compute_30,code=sm_30 -gencode arch=compute_52,code=sm_52 -m64 -odir "src/Device" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-7.5/bin/nvcc -lineinfo -O3 -ftz true -v -Xcompiler -fPIC -Xptxas -dlcm=ca --compile --relocatable-device-code=true -gencode arch=compute_30,code=compute_30 -gencode arch=compute_52,code=compute_52 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_52,code=sm_52 -m64  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


