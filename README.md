# HaploHyped VarAwareML

High-performance genomic data processing pipeline for machine learning. Converts VCF files to optimized HDF5 format with C++ acceleration and Blosc2 compression.

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Key Features

- **559K variants/sec** - C++ VCF parsing with vcfpp
- **6.5x compression** - Blosc2 with LZ4/Zstandard
- **PyTorch ready** - Custom Dataset classes with on-the-fly haplotype encoding
- **GPU support** - CUDA-accelerated processing

## Quick Start

```bash
git clone https://github.com/Jaureguy760/HaploHyped-VarAwareML.git
cd HaploHyped-VarAwareML
./setup.sh
conda activate HaploHyped-VarAwareML
./build.sh
```

## Usage

### VCF to HDF5

```bash
vcf_to_h5 \
    --cohort_name my_study \
    --vcf /path/to/vcf_files \
    --outdir /path/to/output \
    --sample_list samples.txt \
    --cores 10
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

Whole genome (3M variants): ~6s parse, ~12s write

## Architecture

```
VCF Input → C++ Parser (vcfpp) → NumPy Arrays → HDF5 (Blosc2) → PyTorch Dataset
```

## Project Structure

```
├── cpp/           # C++ VCF parser (vcfpp, pybind11)
├── src/
│   ├── haplohyped/   # Core package
│   └── datasets/     # PyTorch datasets
├── tests/         # Test suite
└── docs/          # Documentation
```

## Requirements

- Linux (Ubuntu 20.04+)
- Conda/Mamba
- CUDA (optional)

## Testing

```bash
pytest tests/ -v
```

## License

MIT
