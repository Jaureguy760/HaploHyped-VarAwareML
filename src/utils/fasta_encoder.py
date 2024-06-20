import numpy as np
import h5py
from pysam import FastaFile
from .common_utils import nucleotide_to_index, bitpack_indices, parse_encode_dict

class ReferenceGenome:
    """
    A class to handle encoding and storing reference genome sequences from a FASTA file into HDF5 format.
    
    Attributes:
    fasta_file (str): Path to the input FASTA file.
    encode_spec (dict): Encoding specification dictionary.
    onehot_dict (dict): Dictionary to store encoded sequences for each chromosome.
    """
    
    def __init__(self, fasta_file, encode_spec=None):
        """
        Initialize the ReferenceGenome with the path to the FASTA file and encoding specifications.

        Parameters:
        fasta_file (str): Path to the input FASTA file.
        encode_spec (dict, optional): Encoding specification dictionary. Default is None.
        """
        self.fasta_file = fasta_file
        self.encode_spec = self.parse_encode_dict(encode_spec)
        self.onehot_dict = self.encode_from_fasta(fasta_file)
    
    def encode_from_fasta(self, in_fasta, ignore_case=True):
        """
        Encode sequences from a FASTA file into a dictionary of one-hot encoded sequences.

        Parameters:
        in_fasta (str): Path to the input FASTA file.
        ignore_case (bool, optional): Flag to ignore case when fetching sequences. Default is True.

        Returns:
        dict: A dictionary with chromosome names as keys and encoded sequences as values.
        """
        onehot_dict = {}
        with FastaFile(in_fasta) as fasta:
            for chrom in fasta.references:
                fasta_seq = fasta.fetch(chrom).upper() if ignore_case else fasta.fetch(chrom)
                indices = nucleotide_to_index(fasta_seq)
                packed_indices = bitpack_indices(indices)
                onehot_dict[chrom] = packed_indices
        return onehot_dict
    
    def save_to_h5(self, out_h5):
        """
        Save the encoded sequences to an HDF5 file with compression.

        Parameters:
        out_h5 (str): Path to the output HDF5 file.
        """
        with h5py.File(out_h5, 'w') as f:
            for chrom, encoded_seq in self.onehot_dict.items():
                f.create_dataset(chrom, data=encoded_seq, compression='lzf')
