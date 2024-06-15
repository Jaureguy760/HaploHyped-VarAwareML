import os

# Define the content for each __init__.py file
init_contents = {
    "src": '''"""
HaploHyped-VarAwareML

This package provides tools for processing genomic data, converting VCF files to HDF5 format,
and performing on-the-fly haplotype encoding for machine learning.
"""

# Optionally, you can import important submodules here for easier access
# from .utils import vcf_to_h5, fasta_encoder, h5_reader
# from .datasets import haplotype_dataset
# from .dataloaders import dataloader
# from .models import model
# from .training import train
''',
    "src/utils": '''"""
Utility functions and classes for the HaploHyped-VarAwareML package.
"""

# Import utility modules
from .vcf_to_h5 import read_sample_list, store_individuals, genotype_VCF2hdf5, process_chromosome, merge_h5_files
from .fasta_encoder import ReferenceGenome
from .common_utils import nucleotide_to_index, bitpack_indices
''',
    "src/datasets": '''"""
Dataset classes for the HaploHyped-VarAwareML package.
"""

from .haplotype_dataset import HaplotypeDataset
''',
    "src/dataloaders": '''"""
DataLoader creation functions for the HaploHyped-VarAwareML package.
"""

from .dataloader import get_dataloader
''',
    "src/models": '''"""
Model definitions for the HaploHyped-VarAwareML package.
"""

from .model import YourModelClass
''',
    "src/training": '''"""
Training scripts for the HaploHyped-VarAwareML package.
"""

from .train import train_model
''',
    "tests": '''"""
Unit tests for the HaploHyped-VarAwareML package.
"""
''',
}

# Function to write __init__.py files
def write_init_files(base_dir, init_contents):
    for dir_path, init_content in init_contents.items():
        full_path = os.path.join(base_dir, dir_path)
        init_file_path = os.path.join(full_path, "__init__.py")

        # Create the directory if it doesn't exist
        os.makedirs(full_path, exist_ok=True)

        # Write the content to the __init__.py file
        with open(init_file_path, 'w') as init_file:
            init_file.write(init_content)
        print(f"Updated: {init_file_path}")

# Base directory of your project
base_dir = "/iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML"

# Write the __init__.py files
write_init_files(base_dir, init_contents)
