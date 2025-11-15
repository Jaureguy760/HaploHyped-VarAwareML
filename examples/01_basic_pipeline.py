"""
Basic Pipeline Example - HaploHyped VarAwareML

This example demonstrates the complete workflow from VCF processing to ML training.
"""

import os
from pathlib import Path

# Configuration
DATA_DIR = Path("tests/data")
OUTPUT_DIR = Path("output/example_01")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*70)
print("HaploHyped VarAwareML - Basic Pipeline Example")
print("="*70)
print()

# Step 1: Convert VCF to HDF5
print("Step 1: Converting VCF to HDF5...")
print("-" * 70)

from haplohyped.vcf_to_h5 import VCFtoHDF5Converter

converter = VCFtoHDF5Converter(
    cohort_name="example_cohort",
    vcf_dir=str(DATA_DIR),
    out_dir=str(OUTPUT_DIR),
    sample_list_path=str(DATA_DIR / "ipscs_samples_test.txt"),
    cores=2,
    cxx_threads=1
)

print(f"Processing {len(converter.donor_ids)} samples...")
print(f"Output: {OUTPUT_DIR / 'example_cohort.h5'}")
print()

# Uncomment to run (requires C++ module built):
# converter.run()
print("✓ VCF to HDF5 conversion complete (skipped in example)")
print()

# Step 2: Encode reference genome
print("Step 2: Encoding reference genome...")
print("-" * 70)

from haplohyped.fasta_encoder import ReferenceGenome

ref_output = OUTPUT_DIR / "reference_genome.h5"
print(f"Output: {ref_output}")
print()

# Uncomment to run:
# ref_genome = ReferenceGenome(
#     fasta_file=str(DATA_DIR / "chr22.fasta"),
#     output_dir=str(OUTPUT_DIR / "tmp_ref")
# )
# ref_genome.load_genome_parallel()
print("✓ Reference encoding complete (skipped in example)")
print()

# Step 3: Create PyTorch dataset
print("Step 3: Creating PyTorch dataset...")
print("-" * 70)

from datasets.haplotype_dataset import RandomHaplotypeDataset
from torch.utils.data import DataLoader

# Example configuration (would require HDF5 files from steps 1-2)
print("Dataset configuration:")
print(f"  BED file: {DATA_DIR / 'test_regions.bed'}")
print(f"  Genotype file: {OUTPUT_DIR / 'example_cohort.h5'}")
print(f"  Reference file: {ref_output}")
print(f"  Sequence length: 1000 bp")
print(f"  Batch size: 32")
print()

# Uncomment to create dataset (requires HDF5 files):
# dataset = RandomHaplotypeDataset(
#     bed_file=str(DATA_DIR / "test_regions.bed"),
#     hdf5_genotype_file=str(OUTPUT_DIR / "example_cohort.h5"),
#     hdf5_reference_file=str(ref_output),
#     samples_file=str(DATA_DIR / "ipscs_samples_test.txt"),
#     seq_length=1000,
#     batch_size=32
# )
#
# dataloader = DataLoader(dataset, batch_size=8, num_workers=2)
#
# # Training loop
# for epoch in range(3):
#     for batch_idx, (hap1, hap2) in enumerate(dataloader):
#         print(f"Epoch {epoch}, Batch {batch_idx}: hap1 shape={hap1.shape}, hap2 shape={hap2.shape}")
#         break  # Just show one batch
#     break

print("✓ Dataset setup complete (skipped in example)")
print()

print("="*70)
print("Example Complete!")
print("="*70)
print()
print("To run this example with real data:")
print("1. Build the C++ module: ./build.sh")
print("2. Ensure you have VCF and FASTA files")
print("3. Uncomment the processing steps above")
print("4. Run: python examples/01_basic_pipeline.py")
