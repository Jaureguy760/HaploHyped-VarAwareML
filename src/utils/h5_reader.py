import h5py
import numpy as np

class HDF5Reader:
    def __init__(self, h5_file):
        self.h5_file = h5_file
    
    def fetch_genotypes(self, donor_id, chromosome):
        with h5py.File(self.h5_file, 'r') as f:
            group_path = f'donor_{donor_id}/chr_{chromosome}'
            if group_path in f:
                genotype_data = f[group_path]['genotype'][()]
                return genotype_data
            else:
                raise KeyError(f"No data found for {group_path}")
