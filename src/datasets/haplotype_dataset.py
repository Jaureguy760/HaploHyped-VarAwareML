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
        pos_df = pl.DataFrame(genotype_data, columns=['chrom', 'start', 'stop', 'ref', 'alt', 'phase1', 'phase2'])
        pos_df = pos_df.with_columns(
            pl.col("ref").map_dict(self.encode_spec),
            pl.col("alt").map_dict(self.encode_spec)
        )
        pos_df = pos_df.with_columns(
            pl.when(pl.col("phase1") == 1).then(pl.col("alt")).otherwise(pl.col("ref")).alias("p1_pos"),
            pl.when(pl.col("phase2") == 1).then(pl.col("alt")).otherwise(pl.col("ref")).alias("p2_pos")
        )
        pos_array = pos_df["start"].to_numpy()
        p1_array = pos_df["p1_pos"].to_numpy()
        p2_array = pos_df["p2_pos"].to_numpy()

        hap1 = pl.repeat(ref_sequence, 2).to_numpy()[0]
        hap2 = pl.repeat(ref_sequence, 2).to_numpy()[1]
        hap1[pos_array] = 0
        hap2[pos_array] = 0
        hap1[pos_array, p1_array] = 1
        hap2[pos_array, p2_array] = 1
        return hap1, hap2
