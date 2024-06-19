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



## How to Run the Pipeline

### Data Processing Script

The main script for converting VCF data to HDF5 format and processing the data is `process_vcf.py`. This script can handle storing individuals and chunk-wise processing.

### Usage

1. **Store Individuals:**
   
    This step stores the individual sample names into an HDF5 file with compression.

    ```bash
    python process_vcf.py --cohort_name <cohort_name> --vcf <vcf_directory> --outdir <output_directory> --flag individuals --sample_list <sample_list_file>
    ```

2. **Process VCF to HDF5 (Chunk-wise):**
   
    This step processes the VCF file and stores the genotype data into HDF5 files with compression.

    ```bash
    python process_vcf.py --cohort_name <cohort_name> --vcf <vcf_directory> --outdir <output_directory> --flag chunk --sample_list <sample_list_file>
    ```

### Arguments

- `--cohort_name`: Name of the cohort (required).
- `--vcf`: Path to the directory containing VCF files (required).
- `--outdir`: Path to the directory to save the results (required).
- `--flag`: Type of data to process (`individuals` or `chunk`) (required).
- `--sample_list`: Path to the sample list file (default: `sample_list.txt`).

### Example Commands

**Store Individuals:**

```bash
python process_vcf.py --cohort_name my_study --vcf /path/to/vcf_files --outdir /path/to/output --flag individuals --sample_list sample_list.txt
