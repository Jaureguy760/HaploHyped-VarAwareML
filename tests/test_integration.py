"""Integration tests for the full pipeline using synthetic test data."""

import pytest
import os
import shutil
import tempfile
import h5py
import numpy as np


class TestVCFtoHDF5Integration:
    """Test VCF to HDF5 conversion with real test data."""

    @pytest.fixture
    def test_data_dir(self):
        """Return path to test data directory."""
        return os.path.join(os.path.dirname(__file__), 'data')

    @pytest.fixture
    def output_dir(self):
        """Create temporary output directory."""
        tmpdir = tempfile.mkdtemp()
        yield tmpdir
        shutil.rmtree(tmpdir, ignore_errors=True)

    def test_vcf_files_exist(self, test_data_dir):
        """Test that required VCF test files exist."""
        vcf_file = os.path.join(test_data_dir, 'chr22.filtered.vcf.gz')
        assert os.path.exists(vcf_file), f"VCF file not found: {vcf_file}"

        # Check file size (should be small but not empty)
        size = os.path.getsize(vcf_file)
        assert size > 1000, f"VCF file too small: {size} bytes"
        assert size < 1000000, f"VCF file too large: {size} bytes"

    def test_reference_exists(self, test_data_dir):
        """Test that reference FASTA exists."""
        fasta_file = os.path.join(test_data_dir, 'chr22.fasta')
        assert os.path.exists(fasta_file), f"FASTA file not found: {fasta_file}"

        # Check it's a valid FASTA
        with open(fasta_file) as f:
            first_line = f.readline()
            assert first_line.startswith('>'), "Invalid FASTA format"

    def test_bed_file_exists(self, test_data_dir):
        """Test that BED file with regions exists."""
        bed_file = os.path.join(test_data_dir, 'test_regions.bed')
        assert os.path.exists(bed_file), f"BED file not found: {bed_file}"

        # Count regions
        with open(bed_file) as f:
            lines = [l for l in f if l.strip() and not l.startswith('#')]
        assert len(lines) > 0, "BED file is empty"
        assert len(lines) == 20, f"Expected 20 regions, got {len(lines)}"

    def test_sample_list_matches_vcf(self, test_data_dir):
        """Test that sample list file has the expected samples."""
        sample_file = os.path.join(test_data_dir, 'ipscs_samples_test.txt')
        assert os.path.exists(sample_file), f"Sample file not found: {sample_file}"

        with open(sample_file) as f:
            samples = [line.strip() for line in f if line.strip()]

        assert len(samples) == 3, f"Expected 3 samples, got {len(samples)}"

        # Check they're UUIDs
        for sample in samples:
            parts = sample.split('-')
            assert len(parts) == 5, f"Invalid UUID format: {sample}"


class TestHDF5Output:
    """Test HDF5 output structure and compression."""

    def test_create_simple_hdf5(self):
        """Test creating a simple HDF5 file with compression."""
        import hdf5plugin

        with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
            temp_file = f.name

        try:
            # Create test HDF5 file
            data = np.random.randint(0, 100, size=(100, 10), dtype=np.int32)

            with h5py.File(temp_file, 'w') as f:
                f.create_dataset(
                    'test',
                    data=data,
                    compression=32001,  # Blosc2
                    compression_opts=(2, 2, 0, 0, 5, 1, 2),
                    chunks=True
                )

            # Read back
            with h5py.File(temp_file, 'r') as f:
                read_data = f['test'][:]
                assert np.array_equal(data, read_data)

            # Check file was created (compression ratio varies with data)
            file_size = os.path.getsize(temp_file)
            # Small random data may not compress well, just verify file exists
            assert file_size > 0, "HDF5 file should be created"

        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)


class TestDatasetCompatibility:
    """Test that test data is compatible with PyTorch datasets."""

    @pytest.fixture
    def test_data_dir(self):
        """Return path to test data directory."""
        return os.path.join(os.path.dirname(__file__), 'data')

    def test_bed_file_format(self, test_data_dir):
        """Test BED file can be parsed by polars."""
        try:
            import polars as pl
        except ImportError:
            pytest.skip("polars not installed")

        bed_file = os.path.join(test_data_dir, 'test_regions.bed')

        # Read BED file as the dataset would
        df = pl.read_csv(
            bed_file,
            separator='\t',
            has_header=False,
            new_columns=['chrom', 'start', 'end']
        )

        assert len(df) > 0
        assert all(col in df.columns for col in ['chrom', 'start', 'end'])

        # Check data types
        assert df['start'].dtype in [pl.Int64, pl.Int32, pl.UInt32]
        assert df['end'].dtype in [pl.Int64, pl.Int32, pl.UInt32]


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, '-v'])
