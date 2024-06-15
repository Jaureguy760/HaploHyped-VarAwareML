import torch
from torch.utils.data import Dataset
import numpy as np
import pandas as pd
from .common_utils import unpack_bits, index_to_onehot, parse_encode_dict

class HaplotypeDataset(Dataset):
    def __init__(self, hdf5_reader, donor_ids, chromosomes, reference_genome, encode_spec=None):
        self.hdf5_reader = hdf5_reader
        self.donor_ids = donor_ids
        self.chromosomes = chromosomes
        self.reference_genome = reference_genome
        self.encode_spec = parse_encode_dict(encode_spec)
        self.data = []
        for donor_id in donor_ids:
            for chrom in chromosomes:
                genotype_data = self.hdf5_reader.fetch_genotypes(donor_id, chrom)
                self.data.append((donor_id, chrom, genotype_data))
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        donor_id, chrom, genotype_data = self.data[idx]
        ref_sequence = self.reference_genome.onehot_dict[chrom]
        unpacked_genotypes = unpack_bits(genotype_data)
        hap1, hap2 = self.encode_haplotypes(ref_sequence, unpacked_genotypes)
        return hap1, hap2
    
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

        hap1_onehot = index_to_onehot(hap1, self.encode_spec)
        hap2_onehot = index_to_onehot(hap2, self.encode_spec)

        return torch.tensor(hap1_onehot, dtype=torch.float32), torch.tensor(hap2_onehot, dtype=torch.float32)
