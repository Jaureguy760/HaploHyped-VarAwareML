import numpy as np

# def nucleotide_to_index(seq, encode_spec=None):
#     """
#     Convert a DNA sequence to integer indices.
    
#     Parameters:
#     seq (str): A string representing a DNA sequence.
#     encode_spec (dict, optional): Encoding specification for nucleotides. Defaults to None.
    
#     Returns:
#     np.array: An array of integers representing the indices of nucleotides.
#     """
#     if encode_spec is None:
#         encode_spec = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
#     return np.array([encode_spec.get(nuc, 4) for nuc in seq], dtype=np.int8)

# def bitpack_indices(indices):
#     """
#     Pack nucleotide indices into a 3-bit representation to include `N`.

#     Parameters:
#     indices (np.array): Array of nucleotide indices.

#     Returns:
#     np.array: Packed array of indices in 3-bit representation.
#     """
#     packed = np.packbits(indices.reshape(-1, 2), axis=-1, bitorder='little')
#     return packed


# def index_to_onehot(indices, encode_spec=None):
#     """
#     Convert nucleotide indices to one-hot encoding.
    
#     Parameters:
#     indices (np.array): Array of nucleotide indices.
#     encode_spec (dict, optional): Encoding specification for nucleotides. Defaults to None.
    
#     Returns:
#     np.array: One-hot encoded representation of the indices.
#     """
#     if encode_spec is None:
#         encode_spec = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
#     num_classes = len(encode_spec)
#     return np.eye(num_classes)[indices]
    
# def unpack_bits(packed_data):
#     """
#     Unpack 3-bit packed data back to nucleotide indices.
    
#     Parameters:
#     packed_data (np.array): Packed array of nucleotide indices.
    
#     Returns:
#     np.array: Unpacked array of nucleotide indices.
#     """
#     unpacked = np.unpackbits(packed_data, axis=-1, bitorder='little').reshape(-1, 2)
#     return unpacked


def parse_encode_dict(encode_spec):
    """
    Parse encoding specification into a dictionary.
    
    Parameters:
    encode_spec (str, list, or dict): Encoding specification.
    
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
import numpy as np
import pandas as pd


def array_to_onehot(seq_array, base_list):
    seq_array[np.isin(seq_array, [b"A", b"C", b"G", b"T"], invert=True)] = b"N"  # Convert ambiguous
    return pd.get_dummies(seq_array).reindex(columns=base_list, fill_value=0).to_numpy()

def encode_sequence(seq_data, encode_spec=None, ignore_case=True):
    if isinstance(seq_data, str):
        if ignore_case:
            seq_data = seq_data.upper()
        seq_data = np.fromiter(seq_data, count=len(seq_data), dtype="|S1")
    elif isinstance(seq_data, np.ndarray):
        if seq_data.dtype != "|S1":
            seq_data = seq_data.astype("|S1")
        if ignore_case:
            seq_data = np.char.upper(seq_data)
    else:
        raise TypeError("Please input as string or numpy array!")
    
    encode_spec = parse_encode_dict(encode_spec)
    base_list = list(encode_spec.keys())
    return array_to_onehot(seq_data, base_list)
