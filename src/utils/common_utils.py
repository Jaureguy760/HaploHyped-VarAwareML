import numpy as np

def nucleotide_to_index(seq):
    """
    Convert a DNA sequence to integer indices.
    
    Parameters:
    seq (str): A string representing a DNA sequence.
    
    Returns:
    np.array: An array of integers representing the indices of nucleotides.
    """
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return np.array([mapping[nuc] for nuc in seq], dtype=np.int8)

def bitpack_indices(indices):
    """
    Pack nucleotide indices into a 2-bit representation.
    
    Parameters:
    indices (np.array): Array of nucleotide indices.
    
    Returns:
    np.array: Packed array of indices in 2-bit representation.
    """
    packed = np.packbits(indices.reshape(-1, 4), axis=-1, bitorder='little')
    return packed

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
