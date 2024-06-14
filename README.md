# HaploHyped VarAwareML Pipeline
![HaploHyped VarAware ML PIPELINE Logo](https://github.com/Jaureguy760/HaploHyped-VarAwareML/blob/main/haplohyped.png)


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

## Usage

You can run the pipeline using the command line interface with the following options.

### VCF to HDF5 Conversion

```bash
haplohyped convert --vcf data/example.vcf.gz --cohort_name YourStudy --outdir data/output
```

```
haplohyped train --hdf5 data/output/genotype_data.h5 --bed data/example.bed --epochs 10 --batch_size 32
```

## Model Training

###  VCF to H5
```
haplohyped convert --vcf data/example.vcf.gz --cohort_name YourStudy --outdir data/output

```

### Running Both Stages Sequentially

```
haplohyped convert --vcf data/example.vcf.gz --cohort_name YourStudy --outdir data/output
haplohyped train --hdf5 data/output/genotype_data.h5 --bed data/example.bed --epochs 10 --batch_size 32
```

## Benchmarks and Profiling

### Data Loading Benchmark
```
haplohyped benchmark-data-loading
```

### Encoding Benchmark
```
haplohyped benchmark-encoding
```

### Training Benchmark
```
haplohyped benchmark-training
```

### Profiling
```
haplohyped profile-pipeline
```

## Contributing
```
We welcome contributions to the HaploHyped VarAware ML PIPELINE! See CONTRIBUTING.md for more details.
```

## License
```
This project is licensed under the MIT License - see the LICENSE file for details.
```
## Contact
```
For questions or issues, please open an issue on GitHub.
```
