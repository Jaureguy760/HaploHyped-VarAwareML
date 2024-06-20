import torch
from torch.utils.data import Dataset, DataLoader
import polars as pl
import numpy as np
import random
from mpi4py import MPI
from .common_utils import unpack_bitpacked_data, index_to_onehot, parse_encode_dict
from .h5_reader import HDF5Reader
from .reference_genome import ReferenceGenome

class RandomHaplotypeDataset(Dataset):
    def __init__(self, bed_file, hdf5_reader, reference_genome, encode_spec=None, seed=42, batch_size=1):
        self.bed_df = pl.read_csv(bed_file, sep='\t', has_header=False, new_columns=['chrom', 'start', 'end'])
        self.hdf5_reader = hdf5_reader
        self.reference_genome = reference_genome
        self.encode_spec = parse_encode_dict(encode_spec)
        self.donor_ids = self.hdf5_reader.get_donor_ids()
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.seed = seed
        self.batch_size = batch_size
        self.set_random_seed(self.seed)
        
        # Caching
        self.cached_bed_df = self.bed_df.collect()
        self.cached_reference_genome = self.reference_genome.cache()
        self.cached_hdf5_reader = self.hdf5_reader.cache()
    
    def set_random_seed(self, seed):
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    
    def __len__(self):
        return len(self.cached_bed_df)
    
    def __getitem__(self, idx):
        hap1_batch = []
        hap2_batch = []
        
        for _ in range(self.batch_size):
            while True:
                region_idx = random.randint(0, len(self.cached_bed_df) - 1)
                donor_idx = random.randint(0, len(self.donor_ids) - 1)
                
                region = self.cached_bed_df[region_idx]
                donor_id = self.donor_ids[donor_idx]
                chrom, start, end = region['chrom'], region['start'], region['end']
                
                if not self.cached_reference_genome.has_sequence(chrom, start, end):
                    continue  # Try again with a new random region
                if not self.cached_hdf5_reader.has_genotype(donor_id, chrom, start, end):
                    continue  # Try again with a new random region
                
                ref_sequence = self.cached_reference_genome.get_sequence(chrom, start, end)
                genotype_data = self.cached_hdf5_reader.fetch_genotypes(donor_id, chrom, start, end)
                
                unpacked_genotypes = unpack_bitpacked_data(genotype_data)
                hap1, hap2 = self.encode_haplotypes(ref_sequence, unpacked_genotypes)
                
                hap1_onehot = index_to_onehot(hap1, self.encode_spec)
                hap2_onehot = index_to_onehot(hap2, self.encode_spec)
                
                hap1_batch.append(hap1_onehot)
                hap2_batch.append(hap2_onehot)
                break
        
        hap1_batch = torch.tensor(np.stack(hap1_batch), dtype=torch.float32)
        hap2_batch = torch.tensor(np.stack(hap2_batch), dtype=torch.float32)
        
        return hap1_batch, hap2_batch
    
    def encode_haplotypes(self, ref_sequence, genotype_data):
        # Convert genotype data to NumPy array
        genotype_data = np.array(genotype_data, dtype=np.int8)
        
        # Extract relevant columns from genotype data
        ref = genotype_data[:, 3]
        alt = genotype_data[:, 4]
        phase1 = genotype_data[:, 5]
        phase2 = genotype_data[:, 6]
        
        # Map ref and alt to encode_spec using vectorized operations
        ref_encoded = np.vectorize(self.encode_spec.get)(ref)
        alt_encoded = np.vectorize(self.encode_spec.get)(alt)
        
        # Determine p1_pos and p2_pos using vectorized operations
        p1_pos = np.where(phase1 == 1, alt_encoded, ref_encoded)
        p2_pos = np.where(phase2 == 1, alt_encoded, ref_encoded)
        
        # Get the start positions
        pos_array = genotype_data[:, 1]
        
        # Create copies of ref_sequence using NumPy broadcasting
        hap1 = np.tile(ref_sequence, (len(genotype_data), 1))
        hap2 = np.copy(hap1)
        
        # Update hap1 and hap2 using vectorized operations
        hap1[np.arange(len(genotype_data)), pos_array] = 0
        hap2[np.arange(len(genotype_data)), pos_array] = 0
        hap1[np.arange(len(genotype_data)), pos_array] = p1_pos
        hap2[np.arange(len(genotype_data)), pos_array] = p2_pos
        
        return hap1, hap2
