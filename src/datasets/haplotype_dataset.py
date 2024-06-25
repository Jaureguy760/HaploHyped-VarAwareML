import torch
from torch.utils.data import Dataset, DataLoader
import polars as pl
import numpy as np
import random
import h5py
from .common_utils import unpack_bitpacked_data, index_to_onehot, parse_encode_dict
from .h5_reader import VCFH5Reader
from .reference_genome import ReferenceGenome

def calculate_midpoint_region(start, end, seq_length):
    midpt = (start + end) // 2
    half_seq_length = seq_length // 2
    new_start = max(0, midpt - half_seq_length)
    new_end = midpt + half_seq_length
    return new_start, new_end

class ReferenceGenome:
    def __init__(self, h5_file, encode_spec=None):
        self.h5_file = h5py.File(h5_file, 'r')
        self.encode_spec = parse_encode_dict(encode_spec)

    def get_sequence(self, chrom, start, end):
        sequence = self.h5_file[chrom][start:end]
        return np.array(sequence, dtype="|S1")

    def close(self):
        self.h5_file.close()

class RandomHaplotypeDataset(Dataset):
    def __init__(self, bed_file, hdf5_genotype_file, hdf5_reference_file, samples_file, encode_spec=None, seed=42, batch_size=1, seq_length=1000):
        self.bed_df = pl.read_csv(bed_file, sep='\t', has_header=False, new_columns=['chrom', 'start', 'end'])
        self.vcf_reader = VCFH5Reader(hdf5_genotype_file)
        self.reference_genome = ReferenceGenome(hdf5_reference_file, encode_spec)
        self.encode_spec = parse_encode_dict(encode_spec)
        self.donor_ids = self.read_samples(samples_file)
        self.chromosomes = np.arange(1, 23)
        self.batch_size = batch_size
        self.seq_length = seq_length
        self.set_random_seed(seed)
        self.num_samples = len(self.bed_df)

    def read_samples(self, samples_file):
        with open(samples_file, 'r') as f:
            donors = [line.strip() for line in f]
        return donors

    def set_random_seed(self, seed):
        np.random.seed(seed)

    def __len__(self):
        return self.num_samples

    def __getitem__(self, idx):
        hap1_batch = []
        hap2_batch = []

        for _ in range(self.batch_size):
            region_idx = np.random.randint(0, self.num_samples)
            donor_idx = np.random.randint(0, len(self.donor_ids))
            chrom_idx = np.random.randint(0, len(self.chromosomes))

            region = self.bed_df[region_idx]
            donor_id = self.donor_ids[donor_idx]
            chrom = self.chromosomes[chrom_idx]
            start, end = region['start'], region['end']

            new_start, new_end = calculate_midpoint_region(start, end, self.seq_length)

            ref_sequence = self.reference_genome.get_sequence(chrom, new_start, new_end)
            genotype_data = self.vcf_reader.fetch_genotypes(donor_id, chrom)

            hap1, hap2 = self.encode_haplotypes(ref_sequence, genotype_data, new_start, new_end)

            hap1_onehot = encode_sequence(hap1, self.encode_spec)
            hap2_onehot = encode_sequence(hap2, self.encode_spec)

            hap1_batch.append(hap1_onehot)
            hap2_batch.append(hap2_onehot)

        hap1_batch = np.stack(hap1_batch)
        hap2_batch = np.stack(hap2_batch)

        return torch.tensor(hap1_batch, dtype=torch.float32), torch.tensor(hap2_batch, dtype=torch.float32)

    def encode_haplotypes(self, ref_sequence, genotype_data, start, end):
        seq_length = end - start
        ref_sequence = np.asarray(ref_sequence, dtype="|S1")
        genotype_data = genotype_data.to_numpy(dtype="|S1")

        ref = genotype_data[:, 3]
        alt = genotype_data[:, 4]
        phase1 = genotype_data[:, 5].astype(np.int8)
        phase2 = genotype_data[:, 6].astype(np.int8)

        ref_encoded = np.vectorize(self.encode_spec.get)(ref)
        alt_encoded = np.vectorize(self.encode_spec.get)(alt)

        p1_pos = np.where(phase1 == 1, alt_encoded, ref_encoded)
        p2_pos = np.where(phase2 == 1, alt_encoded, ref_encoded)

        pos_array = genotype_data[:, 1].astype(np.int32) - start

        hap1 = np.zeros((seq_length,), dtype=np.int8)
        hap2 = np.zeros((seq_length,), dtype=np.int8)

        np.put_along_axis(hap1, pos_array, p1_pos, axis=0)
        np.put_along_axis(hap2, pos_array, p2_pos, axis=0)

        return hap1, hap2

    def close(self):
        self.vcf_reader.close()
        self.reference_genome.close()

# Define your parameters
bed_file = 'path/to/bed/file'
vcf_h5_file = 'path/to/vcf/h5file'
reference_h5_file = 'path/to/reference/h5file'
samples_file = 'path/to/samples.txt'
encode_spec = 'path/to/encode/spec'
seed = 42
batch_size = 1
num_epochs = 10
num_workers = 4
seq_length = 1000

# Create the dataset and dataloader
dataset = RandomHaplotypeDataset(bed_file, vcf_h5_file, reference_h5_file, samples_file, encode_spec, seed, batch_size, seq_length)

try:
    dataloader = DataLoader(dataset, batch_size=batch_size, num_workers=num_workers)

    # Train the model using the dataloader
    for epoch in range(num_epochs):
        for batch in dataloader:
            # Model training logic here
            pass

finally:
    # Ensure the files are closed when done
    dataset.close()
