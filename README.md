# HaploHyped VarAwareML Pipeline

<div align="center">

<img src="haplohyped.png" alt="HaploHyped VarAware ML Pipeline Logo" width="500"/>

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)](#testing)
[![HDF5 Compression](https://img.shields.io/badge/compression-Blosc2%20%7C%206.5x-orange.svg)](#performance)

**High-Performance Genomic Data Processing Pipeline for Machine Learning**

[Features](#features) • [Installation](#installation) • [Usage](#usage) • [Performance](#performance) • [Documentation](#documentation) • [Contributing](#contributing)

</div>

---

## 🌟 Overview

HaploHyped VarAwareML is a production-ready, end-to-end pipeline for processing genomic variant data (VCF files) into optimized HDF5 format for machine learning applications. Built with performance and scalability in mind, it leverages C++ for VCF parsing and Blosc2 compression for efficient storage and retrieval.

**Developed by:** Jeff Jaureguy

### Key Highlights

- 🚀 **Blazing Fast**: Processes 500K+ variants/second
- 💾 **Efficient Storage**: 6.5x compression ratio with Blosc2
- 🧬 **ML-Ready**: Direct PyTorch dataset integration
- ⚡ **GPU Accelerated**: On-the-fly haplotype encoding
- 🔧 **Production-Grade**: Comprehensive testing & validation
- 📊 **Scalable**: Handles whole-genome datasets efficiently

---

## 🎯 Features

### Core Functionality

- **VCF to HDF5 Conversion**
  - Parallel processing with configurable threading
  - C++-accelerated VCF parsing using `vcfpp`
  - Automatic chunking and compression
  - Preserves phasing information

- **Reference Genome Encoding**
  - One-hot encoding of reference genomes
  - Parallel chromosome processing
  - Efficient storage with Blosc2 compression

- **PyTorch Integration**
  - Custom `Dataset` classes for on-the-fly haplotype generation
  - Automatic variant application to reference sequences
  - Random sampling for training
  - Multi-worker data loading support

### Technical Features

- **High-Performance Compression**: Blosc2 with LZ4/Zstandard algorithms
- **Memory Efficient**: Streaming processing, chunked I/O
- **Type Safety**: Structured NumPy arrays for genomic data
- **Logging**: Comprehensive logging for production debugging
- **Testing**: Full test suite with 100% core functionality coverage

---

## 📦 Installation

### Prerequisites

- **Conda** or **Mamba** (recommended for faster installation)
- **CUDA** (optional, for GPU acceleration)
- **Linux** (tested on Ubuntu 20.04+)

### Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/Jaureguy760/HaploHyped-VarAwareML.git
cd HaploHyped-VarAwareML

# 2. Set up conda environment
./setup.sh

# 3. Activate environment and build C++ module
conda activate HaploHyped-VarAwareML
./build.sh

# 4. Verify installation
./run_tests.sh
```

### Manual Installation

```bash
# Create environment
conda env create -f environment.yml
conda activate HaploHyped-VarAwareML

# Install Python package
pip install -e .

# Build C++ module
cd cpp && mkdir build && cd build
cmake .. && make
cd ../..
```

---

## 🚀 Usage

### Command Line Interface

#### VCF to HDF5 Conversion

```bash
vcf_to_h5 \
    --cohort_name my_study \
    --vcf /path/to/vcf_files \
    --outdir /path/to/output \
    --sample_list samples.txt \
    --cores 10 \
    --cxx_threads 4
```

**Arguments:**
- `--cohort_name`: Name for your cohort/study
- `--vcf`: Directory containing VCF files (expects `chr{1-22}.filtered.vcf.gz`)
- `--outdir`: Output directory for HDF5 files
- `--sample_list`: Text file with sample IDs (one per line)
- `--cores`: Number of CPU cores for parallel processing
- `--cxx_threads`: Threads for C++ VCF parsing

#### Reference Genome Encoding

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

# Create dataset
dataset = RandomHaplotypeDataset(
    bed_file='regions.bed',
    hdf5_genotype_file='cohort.h5',
    hdf5_reference_file='reference_genome.h5',
    samples_file='samples.txt',
    seq_length=1000,
    batch_size=32
)

# Create dataloader
dataloader = DataLoader(dataset, batch_size=8, num_workers=4)

# Training loop
for epoch in range(num_epochs):
    for hap1, hap2 in dataloader:
        # Your model training here
        predictions = model(hap1, hap2)
        loss = criterion(predictions, targets)
        loss.backward()
        optimizer.step()
```

---

## ⚡ Performance

### Benchmarks

Measured on test system (Intel Xeon, 32 cores, NVMe SSD):

| Operation | Speed | Details |
|-----------|-------|---------|
| VCF Parsing | 559K variants/sec | C++ with vcfpp |
| HDF5 Write | 256K records/sec | With Blosc2 compression |
| HDF5 Read | 342K records/sec | Faster than writing |
| Compression Ratio | 6.5x | Blosc2 LZ4 (level 2) |
| Random Access | 1,597 slices/sec | Sub-millisecond latency |

### Real-World Performance

- **1M variants**: ~2 seconds to parse, ~4 seconds to write HDF5
- **Whole genome** (3M variants): ~6 seconds to parse, ~12 seconds to write
- **100 samples, 22 chromosomes**: ~40 MB compressed (vs 260 MB raw)

### Optimization Features

- ✅ Multi-threaded VCF parsing
- ✅ Parallel chromosome processing
- ✅ Blosc2 compression with multi-threading
- ✅ Optimized slicing with `b2h5py.auto`
- ✅ Memory-mapped file access
- ✅ Chunked HDF5 storage

---

## 📚 Documentation

### Complete Guides

- **[HDF5 Compression Guide](docs/HDF5_COMPRESSION.md)** - Understanding Blosc2 compression
- **[Test Data Documentation](tests/data/README.md)** - Using synthetic test data
- **[Contributing Guidelines](CONTRIBUTING.md)** - How to contribute
- **[Test Results](TEST_RESULTS.md)** - Validation and benchmarks

### Architecture

```
Input (VCF)
    ↓
C++ Parser (vcfpp) ────→ Parallel Processing
    ↓
Structured Arrays (NumPy)
    ↓
HDF5 Writer (Blosc2 Compression) ────→ Optimized Storage
    ↓
PyTorch Dataset ────→ On-the-fly Haplotype Encoding
    ↓
ML Model Training
```

### Project Structure

```
HaploHyped-VarAwareML/
├── cpp/                    # C++ VCF parser
│   ├── parse_vcf.cpp      # Main implementation
│   ├── vcfpp.h            # vcfpp library
│   └── CMakeLists.txt     # Build configuration
├── src/
│   ├── haplohyped/        # Main package
│   │   ├── vcf_to_h5.py   # VCF→HDF5 converter
│   │   ├── fasta_encoder.py  # Reference encoder
│   │   └── main.py        # CLI interface
│   ├── datasets/          # PyTorch datasets
│   │   └── haplotype_dataset.py
│   └── utils/             # Utility functions
├── tests/                 # Test suite
│   ├── data/             # Synthetic test data
│   ├── test_*.py         # Unit & integration tests
│   └── verify_h5.ipynb   # Validation notebook
├── docs/                  # Documentation
└── setup.py              # Package configuration
```

---

## 🧪 Testing

### Run Tests

```bash
# Quick test
pytest tests/ -v

# Full test suite
./run_tests.sh

# Specific test
pytest tests/test_integration.py::TestVCFtoHDF5Integration -v

# With coverage
pytest tests/ --cov=src --cov-report=html
```

### Test Coverage

- ✅ VCF parsing and validation
- ✅ HDF5 compression and integrity
- ✅ Data type preservation
- ✅ Count field accuracy (AC/AF/AN)
- ✅ PyTorch dataset compatibility
- ✅ Performance benchmarks

All tests passing: **6/6 integration tests**, **100% core functionality**

---

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run formatters
black src/ tests/
flake8 src/ tests/

# Run type checker
mypy src/
```

---

## 📊 Use Cases

- **Variant Effect Prediction**: Train models to predict functional impact of variants
- **Population Genetics**: Analyze haplotype diversity across populations
- **Disease Association**: Identify genomic regions associated with phenotypes
- **Functional Genomics**: Integrate with epigenetic data for regulatory analysis
- **Pharmacogenomics**: Predict drug response from genetic variants

---

## 🏆 Acknowledgments

- **vcfpp**: Fast VCF/BCF file parsing ([GitHub](https://github.com/Zilong-Li/vcfpp))
- **pybind11**: Seamless Python-C++ integration ([GitHub](https://github.com/pybind/pybind11))
- **Blosc2**: High-performance compression ([Website](https://www.blosc.org/))
- **HDF Group**: HDF5 file format and libraries

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📧 Contact & Support

- **Issues**: [GitHub Issues](https://github.com/Jaureguy760/HaploHyped-VarAwareML/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Jaureguy760/HaploHyped-VarAwareML/discussions)

---

## 🌟 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{haplohyped2024,
  title={HaploHyped VarAwareML: High-Performance Genomic Data Pipeline},
  author={Ho, Aaron and Jaureguy, Jeff},
  year={2024},
  url={https://github.com/Jaureguy760/HaploHyped-VarAwareML}
}
```

---

<div align="center">

**Built with ❤️ for the genomics and machine learning community**

[⬆ Back to Top](#haplohyped-varawareml-pipeline)

</div>
