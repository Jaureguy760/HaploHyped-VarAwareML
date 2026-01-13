# HaploHyped VarAwareML

High-performance genomic data processing pipeline for machine learning. Converts VCF files to optimized HDF5 format with C++ acceleration and Blosc2 compression.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Features

- **559K variants/sec** - C++ VCF parsing with vcfpp
- **6.5x compression** - Blosc2 with LZ4/Zstandard
- **PyTorch integration** - Custom Dataset classes with on-the-fly haplotype encoding
- **GPU support** - CUDA-accelerated processing

## Installation

### Prerequisites

- Linux (tested on Ubuntu 20.04+)
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
- GCC with C++11 support

### Setup

```bash
# Clone repository
git clone https://github.com/Jaureguy760/HaploHyped-VarAwareML.git
cd HaploHyped-VarAwareML

# Create conda environment
conda env create -f environment.yml
conda activate HaploHyped-VarAwareML

# Build C++ components and install
chmod +x build.sh
./build.sh

# Verify installation
pytest tests/ -v
```

### Manual Installation

If you prefer manual setup:

```bash
# Create environment
conda env create -f environment.yml
conda activate HaploHyped-VarAwareML

# Build C++ VCF parser
cd cpp
mkdir -p build && cd build
cmake .. -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
make -j$(nproc)
cd ../..

# Build Python bindings
cd cpp
g++ -O3 -Wall -shared -fPIC -std=c++11 \
    $(python3 -m pybind11 --includes) parse_vcf.cpp \
    -o parse_vcf$(python3-config --extension-suffix) \
    -I"${CONDA_PREFIX}/include" \
    -L"${CONDA_PREFIX}/lib" \
    -lhts -lz -Wl,-rpath,"${CONDA_PREFIX}/lib"
cd ..

# Install package
pip install -e .
```

## Usage

### VCF to HDF5 Conversion

```bash
vcf_to_h5 \
    --cohort_name my_study \
    --vcf /path/to/vcf_files \
    --outdir /path/to/output \
    --sample_list samples.txt \
    --cores 10
```

### Reference Genome Encoding

```bash
fasta_encoder \
    --fasta reference.fasta \
    --outdir /path/to/output \
    --cores 22
```

### Python API

```python
from datasets import RandomHaplotypeDataset
from torch.utils.data import DataLoader

dataset = RandomHaplotypeDataset(
    bed_file='regions.bed',
    hdf5_genotype_file='cohort.h5',
    hdf5_reference_file='reference_genome.h5',
    samples_file='samples.txt',
    seq_length=1000
)

dataloader = DataLoader(dataset, batch_size=32, num_workers=4)

for hap1, hap2 in dataloader:
    predictions = model(hap1, hap2)
```

## Performance

| Operation | Speed |
|-----------|-------|
| VCF Parsing | 559K variants/sec |
| HDF5 Write | 256K records/sec |
| HDF5 Read | 342K records/sec |
| Compression | 6.5x ratio |

**Benchmarks:** Whole genome (3M variants) processes in ~6s parse, ~12s write on Intel Xeon with NVMe SSD.

## Architecture

```
VCF Input → C++ Parser (vcfpp) → NumPy Arrays → HDF5 (Blosc2) → PyTorch Dataset
```

## Project Structure

```
├── cpp/              # C++ VCF parser (vcfpp, pybind11)
│   ├── parse_vcf.cpp
│   ├── vcfpp.h
│   └── CMakeLists.txt
├── src/
│   ├── haplohyped/   # Core package
│   └── datasets/     # PyTorch datasets
├── tests/            # Test suite
├── docs/             # Documentation
└── environment.yml   # Conda dependencies
```

## Testing

```bash
# Run all tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=src --cov-report=html

# Integration tests only
pytest tests/ -m integration -v
```

## Documentation

- [Architecture Overview](docs/ARCHITECTURE.md)

## License

MIT - see [LICENSE](LICENSE)
