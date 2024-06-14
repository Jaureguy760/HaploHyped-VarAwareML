import numpy as np
import h5py
from pysam import FastaFile

class ReferenceGenome:
    def __init__(self, fasta_file, encode_spec=None):
        """
        Initialize the ReferenceGenome object.

        Parameters:
        fasta_file (str): Path to the reference FASTA file.
        encode_spec (dict, list, tuple, or str, optional): Encoding specification for nucleotides. Defaults to None.
        """
        self.fasta_file = fasta_file
        self.encode_spec = self.parse_encode_dict(encode_spec)
        self.onehot_dict = self.encode_from_fasta(fasta_file)
    
    def parse_encode_dict(self, encode_spec):
        """
        Parse the encoding specification.

        Parameters:
        encode_spec (dict, list, tuple, or str, optional): Encoding specification for nucleotides. Defaults to None.

        Returns:
        dict: Parsed encoding specification.
        """
        if not encode_spec:
            return {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
        elif isinstance(encode_spec, (list, tuple, str)):
            return {base: i for i, base in enumerate(encode_spec)}
        elif isinstance(encode_spec, dict):
            return encode_spec
        else:
            raise TypeError("Please input as dict, list or string!")
    
    def nucleotide_to_index(self, seq):
        """
        Convert a DNA sequence to integer indices.

        Parameters:
        seq (str): A string representing a DNA sequence.

        Returns:
        np.array: An array of integers representing the indices of nucleotides.
        """
        mapping = self.encode_spec
        return np.array([mapping.get(nuc, 4) for nuc in seq], dtype=np.int8)
    
    def bitpack_indices(self, indices):
        """
        Pack nucleotide indices into a 2-bit representation.

        Parameters:
        indices (np.array): Array of nucleotide indices.

        Returns:
        np.array: Packed array of indices in 2-bit representation.
        """
        packed = np.packbits(indices.reshape(-1, 4), axis=-1, bitorder='little')
        return packed

    def encode_from_fasta(self, in_fasta, ignore_case=True):
        """
        Encode sequences from a FASTA file into a one-hot dictionary.

        Parameters:
        in_fasta (str): Path to the input FASTA file.
        ignore_case (bool, optional): Convert lowercase bases to upper, defaults to True.

        Returns:
        dict: Dictionary with chromosome keys and encoded sequences as values.
        """
        onehot_dict = {}
        with FastaFile(in_fasta) as fasta:
            for chrom in fasta.references:
                fasta_seq = fasta.fetch(chrom).upper() if ignore_case else fasta.fetch(chrom)
                encoded_seq = self.nucleotide_to_index(fasta_seq)
                packed_seq = self.bitpack_indices(encoded_seq)
                onehot_dict[chrom] = packed_seq
        return onehot_dict

    def save_to_h5(self, out_h5):
        """
        Save the encoded sequences to an HDF5 file.

        Parameters:
        out_h5 (str): Path to the output HDF5 file.
        """
        with h5py.File(out_h5, 'w') as f:
            for chrom, encoded_seq in self.onehot_dict.items():
                f.create_dataset(chrom, data=encoded_seq, compression='lzf')

# Usage example:
# fasta_file = 'path/to/reference_genome.fasta'
# out_h5 = 'path/to/output_file.h5'
# ref_genome = ReferenceGenome(fasta_file)
# ref_genome.save_to_h5(out_h5)
