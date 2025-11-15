"""Test HDF5 Blosc2 compression functionality."""

import pytest
import tempfile
import os
import numpy as np

# Test different import orders and configurations
def test_blosc2_compression_available():
    """Test that Blosc2 compression filter is available."""
    try:
        import hdf5plugin
        import h5py

        # Check if Blosc2 filter is registered
        filters = h5py.h5z.get_filter_info(32001)
        assert filters is not None, "Blosc2 filter (32001) not registered"
        print("✓ Blosc2 filter is registered")
    except ImportError as e:
        pytest.skip(f"Required libraries not installed: {e}")
    except Exception as e:
        pytest.fail(f"Blosc2 filter check failed: {e}")


def test_blosc2_write_read():
    """Test writing and reading with Blosc2 compression."""
    try:
        import hdf5plugin
        import h5py
    except ImportError:
        pytest.skip("hdf5plugin or h5py not installed")

    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
        temp_file = f.name

    try:
        # Create test data
        data = np.random.randint(0, 100, size=(1000, 100), dtype=np.int32)

        # Write with Blosc2 compression
        with h5py.File(temp_file, 'w') as f:
            f.create_dataset(
                'test_data',
                data=data,
                compression=32001,  # Blosc2
                compression_opts=(2, 2, 0, 0, 5, 1, 2),  # (compcode, clevel, shuffle, ...)
                chunks=True
            )

        print(f"✓ Data written with Blosc2 compression")

        # Read back and verify
        with h5py.File(temp_file, 'r') as f:
            read_data = f['test_data'][:]
            assert np.array_equal(data, read_data), "Data mismatch after compression"

            # Check compression filter
            dset = f['test_data']
            compression_info = dset.compression
            assert compression_info == 32001, f"Expected Blosc2 (32001), got {compression_info}"

        print(f"✓ Data read back successfully, compression verified")

        # Check file size (should be smaller than uncompressed)
        file_size = os.path.getsize(temp_file)
        uncompressed_size = data.nbytes
        compression_ratio = uncompressed_size / file_size

        print(f"✓ Compression ratio: {compression_ratio:.2f}x")
        assert compression_ratio > 1.0, "Data should be compressed"

    finally:
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_blosc2_parameters():
    """Test different Blosc2 compression parameters."""
    try:
        import hdf5plugin
        import h5py
    except ImportError:
        pytest.skip("hdf5plugin or h5py not installed")

    # Test the parameters used in vcf_to_h5.py
    compression_opts = (2, 2, 0, 0, 5, 1, 2)
    # Format: (compcode, clevel, shuffle, unused, blocksize, unused, unused)
    # compcode=2: LZ4
    # clevel=2: compression level

    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as f:
        temp_file = f.name

    try:
        data = np.random.randint(0, 255, size=(1000,), dtype=np.uint8)

        with h5py.File(temp_file, 'w') as f:
            f.create_dataset(
                'test',
                data=data,
                compression=32001,
                compression_opts=compression_opts,
                chunks=True
            )

        with h5py.File(temp_file, 'r') as f:
            read_data = f['test'][:]
            assert np.array_equal(data, read_data)

        print(f"✓ Blosc2 compression parameters work correctly")

    finally:
        if os.path.exists(temp_file):
            os.unlink(temp_file)


def test_b2h5py_import():
    """Test if b2h5py.auto can be imported (optional optimization)."""
    try:
        import b2h5py.auto
        print("✓ b2h5py.auto imported successfully")
    except ImportError:
        pytest.skip("b2h5py not installed (optional)")


if __name__ == "__main__":
    print("Testing HDF5 Blosc2 compression...\n")

    print("1. Testing Blosc2 filter availability...")
    try:
        test_blosc2_compression_available()
    except Exception as e:
        print(f"✗ Failed: {e}\n")

    print("\n2. Testing Blosc2 write/read...")
    try:
        test_blosc2_write_read()
    except Exception as e:
        print(f"✗ Failed: {e}\n")

    print("\n3. Testing Blosc2 parameters...")
    try:
        test_blosc2_parameters()
    except Exception as e:
        print(f"✗ Failed: {e}\n")

    print("\n4. Testing b2h5py.auto import...")
    try:
        test_b2h5py_import()
    except Exception as e:
        print(f"✗ Failed: {e}\n")

    print("\n✓ All compression tests completed!")
