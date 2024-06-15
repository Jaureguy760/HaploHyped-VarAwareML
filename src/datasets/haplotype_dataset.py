import torch
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import numpy as np
import random
from mpi4py import MPI
from .common_utils import unpack_bitpacked_data, index_to_onehot, parse_encode_dict
from .hdf5_reader import HDF5Reader
from .reference_genome import ReferenceGenome

class RandomHaplotypeDataset(Dataset):
    def __init__(self, bed_file, hdf5_reader, reference_genome, encode_spec=None):
        self.bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
        self.hdf5_reader = hdf5_reader
        self.reference_genome = reference_genome
        self.encode_spec = parse_encode_dict(encode_spec)
        self.donor_ids = self.hdf5_reader.get_donor_ids()
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
    
    def __len__(self):
        return len(self.bed_df)
    
    def __getitem__(self, idx):
        region_idx = random.randint(0, len(self.bed_df) - 1)
        donor_idx = random.randint(0, len(self.donor_ids) - 1)
        
        region = self.bed_df.iloc[region_idx]
        donor_id = self.donor_ids[donor_idx]
        chrom, start, end = region['chrom'], region['start'], region['end']
        
        if not self.reference_genome.has_sequence(chrom, start, end):
            return self.__getitem__(idx)  # Try again with a new random region

        if not self.hdf5_reader.has_genotype(donor_id, chrom, start, end):
            return self.__getitem__(idx)  # Try again with a new random region
        
        ref_sequence = self.reference_genome.get_sequence(chrom, start, end)
        genotype_data = self.hdf5_reader.fetch_genotypes(donor_id, chrom, start, end)
        
        unpacked_genotypes = unpack_bitpacked_data(genotype_data)
        hap1, hap2 = self.encode_haplotypes(ref_sequence, unpacked_genotypes)
        
        hap1_onehot = index_to_onehot(hap1, self.encode_spec)
        hap2_onehot = index_to_onehot(hap2, self.encode_spec)
        
        return torch.tensor(hap1_onehot, dtype=torch.float32), torch.tensor(hap2_onehot, dtype=torch.float32)
    
    def encode_haplotypes(self, ref_sequence, genotype_data):
        pos_df = pd.DataFrame(genotype_data, columns=['chrom', 'start', 'stop', 'ref', 'alt', 'phase1', 'phase2'])
        pos_df.replace({"ref": self.encode_spec, "alt": self.encode_spec}, inplace=True)

        pos_df["p1_pos"] = np.where(pos_df["phase1"] == 1, pos_df["alt"], pos_df["ref"])
        pos_df["p2_pos"] = np.where(pos_df["phase2"] == 1, pos_df["alt"], pos_df["ref"])

        pos_array = pos_df["start"].to_numpy()
        p1_array = pos_df["p1_pos"].to_numpy()
        p2_array = pos_df["p2_pos"].to_numpy()

        hap1 = np.copy(ref_sequence)
        hap2 = np.copy(ref_sequence)

        hap1[pos_array] = 0
        hap2[pos_array] = 0

        hap1[pos_array, p1_array] = 1
        hap2[pos_array, p2_array] = 1

        return hap1, hap2

    def has_sequence(self, chrom, start, end):
        """
        Check if the reference genome has the specified sequence.
        
        Parameters:
        chrom (str): Chromosome.
        start (int): Start position.
        end (int): End position.
        
        Returns:
        bool: True if the sequence exists, False otherwise.
        """
        return chrom in self.onehot_dict and end <= len(self.onehot_dict[chrom])

    def has_genotype(self, donor_id, chrom, start, end):
        """
        Check if the genotype data exists for the specified region.
        
        Parameters:
        donor_id (str): Donor ID.
        chrom (str): Chromosome.
        start (int): Start position.
        end (int): End position.
        
        Returns:
        bool: True if the genotype data exists, False otherwise.
        """
        try:
            genotype_data = self.hdf5_reader.fetch_genotypes(donor_id, chrom, start, end)
            return True
        except KeyError:
            return False

# Usage
bed_file = 'path/to/bed_file.bed'
hdf5_reader = HDF5Reader('path/to/hdf5')
reference_genome = ReferenceGenome('path/to/reference.fasta')

dataset = RandomHaplotypeDataset(bed_file, hdf5_reader, reference_genome)
dataloader = DataLoader(dataset, batch_size=32, shuffle=True, num_workers=4)

# Training loop
for batch in dataloader:
    hap1_batch, hap2_batch = batch
    # Train your model here
