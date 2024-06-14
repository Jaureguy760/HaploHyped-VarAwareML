import numpy as np
import h5py
from pysam import FastaFile

class ReferenceGenome:
    def __init__(self, fasta_file, encode_spec=None):
        self.fasta_file = fasta_file
        self.encode_spec = self.parse_encode_dict(encode_spec)
        self.onehot_dict = self.encode_from_fasta(fasta_file)
    
    def parse_encode_dict(self, encode_spec):
        if not encode_spec:
            return {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
        elif isinstance(encode_spec, (list, tuple, str)):
            return {base: i for i, base in enumerate(encode_spec)}
        elif isinstance(encode_spec, dict):
            return encode_spec
        else:
            raise TypeError("Please input as dict, list or string!")
    
    def encode_from_fasta(self, in_fasta, ignore_case=True):
        onehot_dict = {}
        with FastaFile(in_fasta) as fasta:
            for chrom in fasta.references:
                fasta_seq = fasta.fetch(chrom).upper() if ignore_case else fasta.fetch(chrom)
                onehot_dict[chrom] = self.encode_sequence(fasta_seq)
        return onehot_dict
    
    def encode_sequence(self, sequence):
        seq_len = len(sequence)
        encoded_seq = np.zeros((seq_len, len(self.encode_spec)), dtype=np.uint8)
        for i, base in enumerate(sequence):
            if base in self.encode_spec:
                encoded_seq[i, self.encode_spec[base]] = 1
        return encoded_seq
    
    def save_to_h5(self, out_h5):
        with h5py.File(out_h5, 'w') as f:
            for chrom, encoded_seq in self.onehot_dict.items():
                f.create_dataset(chrom, data=encoded_seq, compression='gzip', compression_opts=5)
