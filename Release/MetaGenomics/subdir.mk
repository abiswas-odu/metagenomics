################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../MetaGenomics/Dataset.cpp \
../MetaGenomics/Edge.cpp \
../MetaGenomics/HashTable.cpp \
../MetaGenomics/OverlapGraph.cpp \
../MetaGenomics/Read.cpp \
../MetaGenomics/dna.cpp \
../MetaGenomics/main.cpp 

OBJS += \
./MetaGenomics/Dataset.o \
./MetaGenomics/Edge.o \
./MetaGenomics/HashTable.o \
./MetaGenomics/OverlapGraph.o \
./MetaGenomics/Read.o \
./MetaGenomics/dna.o \
./MetaGenomics/main.o 

CPP_DEPS += \
./MetaGenomics/Dataset.d \
./MetaGenomics/Edge.d \
./MetaGenomics/HashTable.d \
./MetaGenomics/OverlapGraph.d \
./MetaGenomics/Read.d \
./MetaGenomics/dna.d \
./MetaGenomics/main.d 


# Each subdirectory must supply rules for building sources it contributes
MetaGenomics/%.o: ../MetaGenomics/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -g3 -Wall -c -fmessage-length=0 -fopenmp -std=c++11 -O3 -lgomp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


