# Examples

This directory contains example workflows and usage patterns for HaploHyped VarAwareML.

## Available Examples

### 01_basic_pipeline.py
**Complete workflow from VCF to ML training**

Demonstrates:
- VCF to HDF5 conversion
- Reference genome encoding
- PyTorch dataset creation
- Basic training loop setup

Usage:
```bash
python examples/01_basic_pipeline.py
```

## Running Examples

Most examples are designed to work with the synthetic test data in `tests/data/`.

### Prerequisites

1. **Install the package:**
   ```bash
   conda activate HaploHyped-VarAwareML
   pip install -e .
   ```

2. **Build C++ module (for VCF processing):**
   ```bash
   ./build.sh
   ```

3. **Verify test data exists:**
   ```bash
   ls tests/data/
   # Should show: chr22.filtered.vcf.gz, chr22.fasta, test_regions.bed, etc.
   ```

### Example Data

The examples use synthetic genomic data from `tests/data/`:
- **chr22.filtered.vcf.gz**: 1,000 phased SNPs, 3 samples
- **chr22.fasta**: 1Mb reference sequence
- **test_regions.bed**: 20 genomic regions
- **ipscs_samples_test.txt**: 3 sample IDs

## Creating Your Own Examples

When creating new examples:

1. **Use clear documentation**: Add docstrings and comments
2. **Handle errors gracefully**: Check for file existence
3. **Show output**: Print progress and results
4. **Keep it simple**: Focus on one concept per example
5. **Test with synthetic data first**: Use `tests/data/`

## Contributing Examples

We welcome new examples! Please:
1. Follow the existing naming convention (`##_descriptive_name.py`)
2. Add documentation to this README
3. Test with the synthetic data
4. Submit a pull request

## Additional Resources

- **Jupyter Notebooks**: Coming soon
- **Advanced Workflows**: Coming soon
- **Benchmarking Scripts**: Coming soon
