#!/bin/bash
# Build script for HaploHyped-VarAwareML C++ components
# Requires: conda environment activated with htslib, pybind11

set -e

echo "Building HaploHyped-VarAwareML..."

# Check conda environment
if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: No conda environment active. Run: conda activate HaploHyped-VarAwareML"
    exit 1
fi

# Set library paths
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH"

# Build C++ module with CMake
echo "Building C++ VCF parser..."
cd cpp
mkdir -p build
cd build
cmake .. -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
make -j$(nproc)
cd ../..

# Build pybind11 module
echo "Building Python bindings..."
cd cpp
g++ -O3 -Wall -shared -fPIC -std=c++11 \
    $(python3 -m pybind11 --includes) parse_vcf.cpp \
    -o parse_vcf$(python3-config --extension-suffix) \
    -I"${CONDA_PREFIX}/include" \
    -L"${CONDA_PREFIX}/lib" \
    -lhts -lz \
    -Wl,-rpath,"${CONDA_PREFIX}/lib"
cd ..

# Install Python package
echo "Installing Python package..."
pip install -e .

echo ""
echo "Build complete! Test with: pytest tests/ -v"
