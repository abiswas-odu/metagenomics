################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Dataset.cpp \
../Edge.cpp \
../HashTable.cpp \
../OverlapGraph.cpp \
../Read.cpp \
../main.cpp 

OBJS += \
./Dataset.o \
./Edge.o \
./HashTable.o \
./OverlapGraph.o \
./Read.o \
./main.o 

CPP_DEPS += \
./Dataset.d \
./Edge.d \
./HashTable.d \
./OverlapGraph.d \
./Read.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


