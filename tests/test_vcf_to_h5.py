"""Tests for VCF to HDF5 conversion."""

import pytest
import tempfile
import os
from haplohyped.vcf_to_h5 import VCFtoHDF5Converter


class TestVCFtoHDF5Converter:
    """Test VCFtoHDF5Converter class."""

    def test_init(self):
        """Test converter initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_file = os.path.join(tmpdir, 'samples.txt')
            with open(sample_file, 'w') as f:
                f.write('sample1\nsample2\n')

            converter = VCFtoHDF5Converter(
                cohort_name='test_cohort',
                vcf_dir='/path/to/vcf',
                out_dir=tmpdir,
                sample_list_path=sample_file,
                cores=2,
                cxx_threads=1
            )

            assert converter.cohort_name == 'test_cohort'
            assert converter.cores == 2
            assert converter.cxx_threads == 1
            assert len(converter.donor_ids) == 2
            assert converter.donor_ids == ['sample1', 'sample2']
            assert os.path.exists(converter.tmp_dir)

    def test_read_sample_list(self):
        """Test reading sample list from file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write('sample1\nsample2\nsample3\n')
            sample_file = f.name

        try:
            converter = VCFtoHDF5Converter(
                cohort_name='test',
                vcf_dir='/path/to/vcf',
                out_dir='/tmp',
                sample_list_path=sample_file,
                cores=1,
                cxx_threads=1
            )
            samples = converter.donor_ids
            assert len(samples) == 3
            assert samples == ['sample1', 'sample2', 'sample3']
        finally:
            os.unlink(sample_file)

    def test_read_sample_list_file_not_found(self):
        """Test that FileNotFoundError is raised for missing file."""
        with pytest.raises(FileNotFoundError):
            VCFtoHDF5Converter(
                cohort_name='test',
                vcf_dir='/path/to/vcf',
                out_dir='/tmp',
                sample_list_path='/nonexistent/file.txt',
                cores=1,
                cxx_threads=1
            )
