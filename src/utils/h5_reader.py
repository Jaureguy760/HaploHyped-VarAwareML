import h5py
import numpy as np
from .common_utils import unpack_bitpacked_data, nucleotide_to_index, parse_encode_dict

class HDF5Reader:
    def __init__(self, h5_file):
        self.h5_file = h5_file
    
    def fetch_genotypes(self, donor_id, chromosome):
        with h5py.File(self.h5_file, 'r') as f:
            group_path = f'donor_{donor_id}/chr_{chromosome}'
            if group_path in f:
                packed_data = f[group_path]['genotype'][()]
                unpacked_data = unpack_bitpacked_data(packed_data)
                return unpacked_data
            else:
                raise KeyError(f"No data found for {group_path}")
