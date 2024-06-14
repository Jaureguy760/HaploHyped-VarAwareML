import os
import pandas as pd
import numpy as np
import h5py
from multiprocessing import Pool

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
    # Pack indices into 2-bit representation
    packed = np.packbits(indices.reshape(-1, 4), axis=-1, bitorder='little')
    return packed

def probes_VCF2hdf5(data_path, save_path, study_name, chunk_size=100000):
    """
    Process VCF file and store probe data into an HDF5 file with compression.
    
    Parameters:
    data_path (str): Path to the input VCF file.
    save_path (str): Path to save the HDF5 file.
    study_name (str): Name of the study.
    chunk_size (int): Number of rows to process at a time.
    """
    h5_file = os.path.join(save_path, 'probes', study_name + '.h5')
    
    # Check if the HDF5 file exists and is up-to-date
    if os.path.isfile(h5_file) and os.path.getmtime(data_path) <= os.path.getmtime(h5_file):
        print(f"{h5_file} is already up-to-date.")
        return

    if os.path.isfile(h5_file):
        os.remove(h5_file)

    hash_table = {'keys': np.array([], dtype=np.int64), 'allele': np.array([], dtype='S')}
    with pd.read_csv(data_path, sep='\t', chunksize=chunk_size, header=None, index_col=None) as df_chunks:
        with h5py.File(h5_file, 'w') as hf:
            for i, chunk in enumerate(df_chunks):
                chunk.columns = ["CHR", "bp", "ID", 'allele1', 'allele2', 'QUAL', 'FILTER', 'INFO']
                allele1_indices = chunk['allele1'].map({'A': 0, 'C': 1, 'G': 2, 'T': 3}).astype(np.int8)
                allele2_indices = chunk['allele2'].map({'A': 0, 'C': 1, 'G': 2, 'T': 3}).astype(np.int8)
                
                # Bitpack the indices
                packed_allele1 = bitpack_indices(allele1_indices)
                packed_allele2 = bitpack_indices(allele2_indices)
                
                packed_chunk = np.column_stack((chunk[["CHR", "bp", "ID"]], packed_allele1, packed_allele2))

                if 'data' not in hf:
                    maxshape = (None, packed_chunk.shape[1])
                    # Using default chunking provided by HDF5
                    dset = hf.create_dataset('data', data=packed_chunk, maxshape=maxshape, compression='lzf')
                else:
                    hf['data'].resize((hf['data'].shape[0] + packed_chunk.shape[0]), axis=0)
                    hf['data'][-packed_chunk.shape[0]:] = packed_chunk

    pd.DataFrame.from_dict(hash_table).to_csv(os.path.join(save_path, 'probes', study_name + '_hash_table.csv.gz'), index=False, compression='gzip', sep='\t')

def ind_VCF2hdf5(data_path, save_path, study_name):
    """
    Process VCF file and store individual data into an HDF5 file with compression.
    
    Parameters:
    data_path (str): Path to the input VCF file.
    save_path (str): Path to save the HDF5 file.
    study_name (str): Name of the study.
    """
    h5_file = os.path.join(save_path, 'individuals', study_name + '.h5')
    
    # Check if the HDF5 file exists and is up-to-date
    if os.path.isfile(h5_file) and os.path.getmtime(data_path) <= os.path.getmtime(h5_file):
        print(f"{h5_file} is already up-to-date.")
        return
    
    if os.path.isfile(h5_file):
        os.remove(h5_file)

    with open(data_path, 'r') as f:
        individuals = [line.strip() for line in f]
    
    with h5py.File(h5_file, 'w') as hf:
        # Using default chunking provided by HDF5
        hf.create_dataset('individuals', data=np.array(individuals, dtype='S'), compression='lzf')

def genotype_VCF2hdf5(data_path, id, save_path, study_name):
    """
    Process VCF file and store genotype data into an HDF5 file with compression.
    
    Parameters:
    data_path (str): Path to the input VCF file.
    id (int): Identifier for the genotype.
    save_path (str): Path to save the HDF5 file.
    study_name (str): Name of the study.
    """
    h5_file = os.path.join(save_path, 'genotype', str(id) + '_' + study_name + '.h5')
    
    # Check if the HDF5 file exists and is up-to-date
    if os.path.isfile(h5_file) and os.path.getmtime(data_path) <= os.path.getmtime(h5_file):
        print(f"{h5_file} is already up-to-date.")
        return

    df = pd.read_csv(data_path, header=None, index_col=None, sep='\t', dtype=str)
    data = df.to_numpy()
    
    with h5py.File(h5_file, 'w') as hf:
        # Using default chunking provided by HDF5
        hf.create_dataset('genotype', data=data.astype(np.int8), compression='lzf')

    os.remove(data_path)

def process_chromosome(chromosome, vcf_dir, out_dir, study_name):
    """
    Process a single chromosome's VCF file and store the data into HDF5 files.
    
    Parameters:
    chromosome (int): Chromosome number to process.
    vcf_dir (str): Directory containing VCF files.
    out_dir (str): Directory to save the output HDF5 files.
    study_name (str): Name of the study.
    """
    vcf_file = os.path.join(vcf_dir, f'chr{chromosome}.filtered.vcf.gz')
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'genotype'), exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'individuals'), exist.ok=True)
    os.makedirs(os.path.join(out_dir, 'probes'), exist.ok=True)
    os.makedirs(os.path.join(out_dir, 'tmp_files'), exist.ok=True)
    
    probes_VCF2hdf5(vcf_file, out_dir, study_name, chunk_size=100000)  # Adjust chunk size here
    ind_VCF2hdf5(vcf_file, out_dir, study_name)

if __name__ == "__main__":
    vcf_dir = "/path/to/vcf_dir"
    out_dir = "/path/to/output_dir"
    study_name = "your_study_name"
    chromosomes = range(1, 23)

    with Pool(processes=4) as pool:
        pool.starmap(process_chromosome, [(chr, vcf_dir, out_dir, study_name) for chr in chromosomes])
