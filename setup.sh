#!/bin/bash

set -e

LOGFILE="setup_and_run.log"
exec > >(tee -i ${LOGFILE})
exec 2>&1

echo "Starting setup and compilation..."

# Create and activate the Conda environment
echo "Creating and activating the Conda environment..."
mamba env create -f environment.yml
source /iblm/netapp/home/jjaureguy/mambaforge/etc/profile.d/conda.sh
conda activate HaploHyped-VarAwareML

# Set LD_LIBRARY_PATH to include the conda lib directory
echo "Setting LD_LIBRARY_PATH..."
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib:$LD_LIBRARY_PATH
# Build the C++ program
echo "Building the C++ program..."
cd cpp
mkdir -p build
cd build
cmake .. || { echo "cmake failed"; exit 1; }
make || { echo "make failed"; exit 1; }
cd ../../

# Compile the Python module using pybind11
echo "Compiling the Python module using pybind11..."
# g++ -shared -fPIC -fsanitize=address -fno-omit-frame-pointer -O1 -g $(python3 -m pybind11 --includes) parse_vcf.cpp -o parse_vcf$(python3-config --extension-suffix) -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -lpython3.8 -lcrypt -lpthread -ldl -lutil -lm -lhts


# g++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) parse_vcf.cpp -o parse_vcf$(python3-config --extension-suffix) -I/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/include -I/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/include/htslib -I/iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML/cpp -L/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib -lhts -llzma -ldeflate -lbz2 -lz -Wl,-rpath,/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib
# Verify the compiled module
echo "Verifying the compiled module..."
ldd cpp/build/parse_vcf$(python3-config --extension-suffix) || { echo "ldd verification failed"; exit 1; }


g++ -O3 -Wall -shared -fPIC -std=c++11 \
    $(python3 -m pybind11 --includes) parse_vcf.cpp \
    -o parse_vcf$(python3-config --extension-suffix) \
    -I/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/include/python3.8 \
    -I/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/include \
    -L/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib \
    -lpython3.8 -lhts -lcrypt -lz





# Set PYTHONPATH to include the directory with the compiled module
echo "Setting PYTHONPATH..."
export PYTHONPATH=/iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML/cpp/build:$PYTHONPATH

echo "Setup and execution complete. The C++ program and Python module have been compiled successfully, and the Python script has been executed."
