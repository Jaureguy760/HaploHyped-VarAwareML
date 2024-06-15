# HaploHyped VarAwareML Pipeline
<img src="https://github.com/Jaureguy760/HaploHyped-VarAwareML/blob/main/haplohyped.png" alt="HaploHyped VarAware ML PIPELINE Logo" width="400"/>

The HaploHyped VarAwareML Pipeline is an integrated, end-to-end solution for processing genomic data, converting VCF files to HDF5 format, and performing GPU-accelerated on-the-fly haplotype encoding for machine learning.

## Features

- End-to-end workflow from raw VCF processing to ML model training
- Parallel processing and optimization techniques
- Efficient data storage with HDF5
- GPU-accelerated on-the-fly haplotype encoding
- Seamless integration with PyTorch for ML

## Installation

Follow these steps to install the HaploHyped VarAwareML Pipeline.

### Prerequisites

- Conda
- CUDA for GPU acceleration
- Nextflow

### Steps

1. Clone the repository:

    ```bash
    git clone https://github.com/Jaureguy760/HaploHyped-VarAwareML.git
    cd HaploHyped-VarAwareML
    ```

2. Create a Conda environment:

    ```bash
    conda env create -f environment.yml
    conda activate haplohyped-env
    ```

3. Install the package:

    ```bash
    pip install -e .
    ```

4. Install Nextflow:
   
    Follow the instructions on the [Nextflow website](https://www.nextflow.io/) to install Nextflow.
## Data Storage

### Reading Data in Chunks
```
The data is read from VCF files in chunks to manage memory usage efficiently.
Each chunk is processed individually.
```
### Nucleotide Encoding
```
Nucleotides are mapped to integer indices using int8 type. For example:
'A' -> 0
'C' -> 1
'G' -> 2
'T' -> 3
```
### Bitpacking
```
The integer indices are then packed into a more compact form using bitpacking. This reduces the amount of space needed to store the data.
For four nucleotides, each represented by 2 bits, they can be packed into a single byte (8 bits).
```
### HDF5 Storage with Compression
The processed data is stored in HDF5 format using the lzf compression algorithm, which provides a good balance between compression efficiency and read/write speed.
HDF5 automatically handles chunking unless explicitly specified, allowing for efficient storage and retrieval.
```
