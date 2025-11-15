# Test Data

This directory contains synthetic test data for the HaploHyped-VarAwareML pipeline.

## Files

### VCF Data
- **chr22.filtered.vcf.gz** - Synthetic VCF file with 1,000 SNP variants
  - Chromosome: chr22
  - Variants: 1,000 phased SNPs (positions 10M-20M)
  - Samples: 3 (matching `ipscs_samples_test.txt`)
  - Format: VCF 4.2, gzip compressed
  - All genotypes are phased (e.g., 0|1, 1|0, etc.)

### Reference Genome
- **chr22.fasta** - Synthetic reference genome sequence
  - Chromosome: chr22
  - Length: 1,000,000 bp
  - Randomly generated sequence (for testing purposes only)

### Regions
- **test_regions.bed** - BED file with 20 test genomic regions
  - Each region is 1,000 bp
  - Covers positions across chr22 (10M-20M)
  - Used for PyTorch dataset testing

### Samples
- **ipscs_samples_test.txt** - Sample IDs (3 samples)
  - Matches the samples in the VCF file
  - One sample ID per line

## Usage

### Test VCF to HDF5 Conversion

```bash
cd /path/to/HaploHyped-VarAwareML

# Create output directory
mkdir -p tests/output

# Run conversion (requires environment activated and C++ module built)
vcf_to_h5 \
    --cohort_name test_cohort \
    --vcf tests/data \
    --outdir tests/output \
    --sample_list tests/data/ipscs_samples_test.txt \
    --cores 2 \
    --cxx_threads 1
```

### Test FASTA Encoding

```bash
fasta_encoder \
    --fasta tests/data/chr22.fasta \
    --outdir tests/output \
    --cores 2
```

### Test PyTorch Dataset

```python
from datasets import RandomHaplotypeDataset
from torch.utils.data import DataLoader

dataset = RandomHaplotypeDataset(
    bed_file='tests/data/test_regions.bed',
    hdf5_genotype_file='tests/output/test_cohort.h5',
    hdf5_reference_file='tests/output/reference_genome.h5',
    samples_file='tests/data/ipscs_samples_test.txt',
    seq_length=1000,
    batch_size=4
)

dataloader = DataLoader(dataset, batch_size=1, num_workers=2)

for hap1, hap2 in dataloader:
    print(f"Haplotype 1 shape: {hap1.shape}")
    print(f"Haplotype 2 shape: {hap2.shape}")
    break
```

## Data Characteristics

| File | Type | Size | Description |
|------|------|------|-------------|
| chr22.filtered.vcf.gz | VCF | ~21KB | 1,000 phased SNPs, 3 samples |
| chr22.fasta | FASTA | ~989KB | 1Mb reference sequence |
| test_regions.bed | BED | 480B | 20 genomic regions |
| ipscs_samples_test.txt | Text | 110B | 3 sample IDs |

## Notes

- This is **synthetic test data** for pipeline testing only
- Not suitable for real biological analysis
- Data is small enough to run quickly on any system
- All positions and sequences are randomly generated
- The VCF follows standard format with phased genotypes (GT field uses `|`)

## Regenerating Test Data

If you need to regenerate the test data:

```bash
cd tests/data
# The data generation scripts are embedded in the repository setup
# See the commit history for the Python scripts used to generate this data
```

## Real Data

For production use, you'll need:
1. Real VCF files from sequencing (typically very large, use `--vcf` directory)
2. Reference genome FASTA (e.g., hg38.fa from UCSC or Ensembl)
3. BED file with regions of interest
4. Sample list matching your VCF
