import numpy as np
import h5py
from pysam import FastaFile
from .common_utils import nucleotide_to_index, bitpack_indices, parse_encode_dict

class ReferenceGenome:
    def __init__(self, fasta_file, encode_spec=None):
        self.fasta_file = fasta_file
        self.encode_spec = self.parse_encode_dict(encode_spec)
        self.onehot_dict = self.encode_from_fasta(fasta_file)
    
    def encode_from_fasta(self, in_fasta, ignore_case=True):
        onehot_dict = {}
        with FastaFile(in_fasta) as fasta:
            for chrom in fasta.references:
                fasta_seq = fasta.fetch(chrom).upper() if ignore_case else fasta.fetch(chrom)
                indices = self.nucleotide_to_index(fasta_seq)
                packed_indices = self.bitpack_indices(indices)
                onehot_dict[chrom] = packed_indices
        return onehot_dict
    
    def save_to_h5(self, out_h5):
        with h5py.File(out_h5, 'w') as f:
            for chrom, encoded_seq in self.onehot_dict.items():
                f.create_dataset(chrom, data=encoded_seq, compression='lzf')
