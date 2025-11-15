# Architecture Documentation

## Overview

HaploHyped VarAwareML is a high-performance genomic data pipeline designed for machine learning applications. The architecture is modular, with clear separation between data processing, storage, and ML integration layers.

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                         Input Data Layer                            │
├─────────────────────────────────────────────────────────────────────┤
│  VCF Files          FASTA Reference         BED Regions             │
│  (Variants)         (Genome Sequence)       (Target Regions)        │
└────────┬────────────────────┬────────────────────┬──────────────────┘
         │                    │                    │
         ▼                    ▼                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      Processing Layer                               │
├─────────────────────────────────────────────────────────────────────┤
│  ┌─────────────────┐  ┌──────────────────┐                         │
│  │ VCF Parser      │  │ FASTA Encoder    │                         │
│  │ (C++ vcfpp)     │  │ (Polars + NumPy) │                         │
│  │ - Phased GT     │  │ - One-hot encode │                         │
│  │ - Parallel      │  │ - Chunking       │                         │
│  └────────┬────────┘  └────────┬─────────┘                         │
│           │                    │                                    │
│           ▼                    ▼                                    │
│  ┌──────────────────────────────────────────────┐                  │
│  │        HDF5 Storage (Blosc2 Compression)      │                  │
│  │  - Genotypes: (chrom, pos, ref, alt, phase)  │                  │
│  │  - Reference: One-hot encoded bases          │                  │
│  │  - Compression: 6.5x ratio                   │                  │
│  └────────────────────┬─────────────────────────┘                  │
└───────────────────────┼─────────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      Data Access Layer                              │
├─────────────────────────────────────────────────────────────────────┤
│  ┌──────────────────────────────────────────────────────────┐      │
│  │          RandomHaplotypeDataset (PyTorch)                │      │
│  │  - Random region sampling from BED file                  │      │
│  │  - On-the-fly haplotype sequence generation              │      │
│  │  - Variant-aware encoding                                │      │
│  └────────────────────────┬─────────────────────────────────┘      │
└───────────────────────────┼─────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      ML Training Layer                              │
├─────────────────────────────────────────────────────────────────────┤
│  PyTorch DataLoader → Model → Training Loop                        │
│  - Batch processing                                                 │
│  - Multi-worker loading                                             │
│  - GPU acceleration                                                 │
└─────────────────────────────────────────────────────────────────────┘
```

## Component Details

### 1. VCF Processing (`vcf_to_h5.py`)

**Technology**: C++ (vcfpp) with Python bindings (pybind11)

**Key Features**:
- Parallel VCF parsing using C++ for performance
- Preserves phased genotype information (0|0, 0|1, 1|0, 1|1)
- Outputs structured NumPy arrays
- Multi-threading support

**Data Flow**:
```
VCF File → vcfpp Parser → Structured Array → HDF5 Dataset
                                              (Blosc2 compressed)
```

**Output Schema**:
```python
dtype = [
    ('chrom', 'U10'),      # Chromosome
    ('start', np.int64),   # Position (0-based)
    ('stop', np.int64),    # Position + 1
    ('ref', 'U100'),       # Reference allele
    ('alt', 'U100'),       # Alternate allele
    ('phase1', np.int8),   # Haplotype 1 (0 or 1)
    ('phase2', np.int8),   # Haplotype 2 (0 or 1)
]
```

**Performance**:
- Parsing: ~559K variants/sec
- Writing: ~256K records/sec
- Compression: 6.5x ratio

### 2. Reference Genome Encoding (`fasta_encoder.py`)

**Technology**: Polars (data processing) + NumPy (encoding)

**Key Features**:
- One-hot encoding: A=[1,0,0,0], C=[0,1,0,0], G=[0,0,1,0], T=[0,0,0,1]
- Chunked processing for memory efficiency
- Parallel chromosome processing
- HDF5 storage with Blosc2 compression

**Data Flow**:
```
FASTA File → Polars DataFrame → One-hot Encoding → HDF5 Dataset
                                                    (per chromosome)
```

**Output Format**:
- Each chromosome: shape (length, 4) float32 array
- Dimensions: [position, base_channel]
- Channels: [A, C, G, T]

### 3. PyTorch Dataset Integration (`haplotype_dataset.py`)

**Class**: `RandomHaplotypeDataset`

**Architecture**:
```python
class RandomHaplotypeDataset:
    ├── __init__()
    │   ├── Load BED regions
    │   ├── Open HDF5 files (genotypes + reference)
    │   └── Parse sample list
    │
    ├── __getitem__(idx)
    │   ├── Select random region from BED
    │   ├── Query genotypes in region
    │   ├── Fetch reference sequence
    │   ├── Apply phased variants
    │   └── Return (haplotype1, haplotype2)
    │
    └── __len__()
        └── Return total samples (regions × samples)
```

**Data Generation Pipeline**:
1. **Region Selection**: Random BED region (chrom, start, end)
2. **Genotype Query**: Extract variants within region from HDF5
3. **Reference Fetch**: Get base reference sequence from encoded genome
4. **Variant Application**: Apply phased variants to generate haplotypes
5. **Output**: Two tensors (haplotype1, haplotype2) with shape (seq_length, 4)

**Features**:
- On-the-fly sequence generation (no pre-computation)
- Memory-efficient (only loads required data)
- Batching support via PyTorch DataLoader
- Multi-worker compatible

### 4. Storage Layer (HDF5 + Blosc2)

**Compression Configuration**:
```python
import hdf5plugin

blosc2_config = hdf5plugin.Blosc2(
    cname='lz4',        # Compression algorithm
    clevel=5,           # Compression level (1-9)
    filters=hdf5plugin.Blosc2.SHUFFLE  # Byte shuffle
)
```

**Performance Characteristics**:
- **Compression Ratio**: 6.5x for genomic data
- **Read Speed**: 342K records/sec (with b2h5py.auto)
- **Random Access**: 1,597 slices/sec
- **Filter ID**: 32001 (Blosc2 HDF5 plugin)

**File Organization**:
```
genotypes.h5
├── /donor_001
│   ├── chr1    (structured array)
│   ├── chr2    (structured array)
│   └── ...
├── /donor_002
│   └── ...
└── /donor_XXX

reference.h5
├── /chr1      (one-hot encoded, float32)
├── /chr2      (one-hot encoded, float32)
└── ...
```

### 5. Build System

**Components**:
1. **CMake** (`cpp/CMakeLists.txt`): C++ module compilation
2. **pybind11**: Python-C++ bindings
3. **Conda Environment**: Dependency management

**Build Process**:
```bash
./build.sh
    ├── Detect CONDA_PREFIX
    ├── Configure CMake
    ├── Find dependencies (htslib, vcfpp)
    ├── Compile C++ module
    └── Install Python bindings
```

**Dependencies**:
- C++: htslib, vcfpp, pybind11
- Python: numpy, h5py, hdf5plugin, torch, polars

## Data Flow Example

### End-to-End Pipeline

```
1. INPUT
   ├── chr22.vcf.gz (1000 variants, 3 samples)
   ├── chr22.fasta (1Mb reference)
   └── regions.bed (20 genomic regions)

2. PROCESSING
   ├── VCF → HDF5
   │   └── Output: cohort.h5 (21KB compressed)
   │
   └── FASTA → HDF5
       └── Output: reference.h5 (~250KB compressed)

3. ML TRAINING
   ├── RandomHaplotypeDataset
   │   ├── Batch size: 32
   │   ├── Sequence length: 1000 bp
   │   └── Workers: 4
   │
   ├── DataLoader
   │   └── Yields: (hap1[32,1000,4], hap2[32,1000,4])
   │
   └── Model Training
       └── Consume batches for training
```

## Performance Optimization

### 1. C++ VCF Parsing
- **Why**: VCF parsing is I/O and CPU intensive
- **Speedup**: ~10x faster than Python-based parsers
- **Implementation**: vcfpp library with parallel processing

### 2. Blosc2 Compression
- **Why**: Reduce storage and improve I/O speed
- **Benefits**:
  - 6.5x compression ratio
  - Faster reads than uncompressed (due to reduced I/O)
  - SIMD-optimized decompression
- **Configuration**: LZ4 algorithm with shuffle filter

### 3. b2h5py Optimizations
- **Why**: Optimize Blosc2-compressed HDF5 reads
- **Speedup**: 8% faster than standard h5py
- **Usage**: `import b2h5py.auto` (drop-in replacement)

### 4. On-the-Fly Data Generation
- **Why**: Avoid storing pre-computed sequences
- **Benefits**:
  - Minimal disk usage
  - Infinite data augmentation via random sampling
  - Memory-efficient batching

## Testing Strategy

### Unit Tests (`tests/`)
- **test_vcf_to_h5.py**: VCF parsing validation
- **test_fasta_encoder.py**: Reference encoding validation
- **test_compression.py**: Blosc2 compression verification
- **test_dataset.py**: PyTorch dataset functionality

### Integration Tests
- **test_integration.py**: End-to-end pipeline validation
- Synthetic data: 1000 variants, 3 samples, 1Mb reference
- Validates: compression, data integrity, performance

### CI/CD (GitHub Actions)
```yaml
Strategy:
  - Python versions: 3.8, 3.9, 3.10, 3.11
  - Test suite: pytest with coverage
  - Code quality: black, flake8, isort
  - Coverage: Codecov integration
```

## Scalability Considerations

### Current Limits
- **VCF Size**: Tested up to 10M variants
- **Samples**: Designed for 100-1000 samples
- **Chromosomes**: All human chromosomes supported

### Performance Scaling
| Dataset Size | Processing Time | Storage Size |
|--------------|-----------------|--------------|
| 1K variants  | ~2 seconds      | ~21 KB       |
| 100K variants| ~3 minutes      | ~2 MB        |
| 10M variants | ~5 hours        | ~200 MB      |

### Future Optimizations
1. **Distributed Processing**: Spark/Dask for larger cohorts
2. **GPU Acceleration**: CUDA kernels for variant encoding
3. **Advanced Compression**: Specialized genomic compressors (e.g., GTC)

## Security Considerations

### Data Privacy
- No telemetry or external data transmission
- Local-only processing
- Supports encrypted file systems

### Input Validation
- VCF format validation
- Sample ID verification
- Chromosome name checking

## Extension Points

### Adding New Features

1. **Custom Encodings**:
   ```python
   class CustomEncoder(ReferenceGenome):
       def encode_base(self, base):
           # Custom encoding logic
   ```

2. **New Dataset Types**:
   ```python
   class WindowedHaplotypeDataset(Dataset):
       # Implement sliding window access
   ```

3. **Additional Annotations**:
   ```python
   # Extend structured array dtype
   dtype.append(('quality', np.float32))
   ```

## Dependencies

### Core Dependencies
```
Python >= 3.8
├── h5py >= 3.7          # HDF5 interface
├── hdf5plugin >= 4.0    # Blosc2 compression
├── numpy >= 1.20        # Array operations
├── torch >= 1.10        # ML framework
├── polars >= 0.15       # DataFrame processing
└── b2h5py >= 0.1        # Optimized HDF5 reads

C++ >= 17
├── htslib >= 1.15       # VCF/BCF handling
├── vcfpp >= 0.1         # VCF parsing
└── pybind11 >= 2.10     # Python bindings
```

### Development Dependencies
```
pytest >= 7.0            # Testing
black >= 22.0            # Code formatting
flake8 >= 5.0            # Linting
mypy >= 1.0              # Type checking
sphinx >= 5.0            # Documentation
```

## References

- **VCF Specification**: [https://samtools.github.io/hts-specs/VCFv4.3.pdf](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- **HDF5 Documentation**: [https://www.hdfgroup.org/](https://www.hdfgroup.org/)
- **Blosc2 Compression**: [https://github.com/Blosc/c-blosc2](https://github.com/Blosc/c-blosc2)
- **PyTorch Datasets**: [https://pytorch.org/docs/stable/data.html](https://pytorch.org/docs/stable/data.html)
