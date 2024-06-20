import h5py
import numpy as np
from .common_utils import unpack_bitpacked_data, nucleotide_to_index, parse_encode_dict

class HDF5Reader:
    """
    A class to read and fetch genotype data from HDF5 files.

    Attributes:
    h5_file (str): Path to the HDF5 file containing genotype data.
    """
    
    def __init__(self, h5_file):
        """
        Initialize the HDF5Reader with the path to the HDF5 file.

        Parameters:
        h5_file (str): Path to the HDF5 file.
        """
        self.h5_file = h5_file
    
    def fetch_genotypes(self, donor_id, chromosome):
        """
        Fetch and unpack genotype data for a specific donor and chromosome.

        Parameters:
        donor_id (str): Identifier for the donor.
        chromosome (int): Chromosome number.

        Returns:
        np.ndarray: Unpacked genotype data.

        Raises:
        KeyError: If no data is found for the specified donor and chromosome.
        """
        with h5py.File(self.h5_file, 'r') as f:
            group_path = f'donor_{donor_id}/chr_{chromosome}'
            if group_path in f:
                packed_data = f[group_path]['genotype'][()]
                unpacked_data = unpack_bitpacked_data(packed_data)
                return unpacked_data
            else:
                raise KeyError(f"No data found for {group_path}")
