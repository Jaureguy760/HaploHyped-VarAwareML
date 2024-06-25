import os
import argparse
import h5py
import numpy as np
from multiprocessing import Pool, cpu_count
import logging
import polars as pl
from pysam import FastaFile
import shutil
import hdf5plugin
import b2h5py.auto

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("reference_genome.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def array_to_onehot(seq_array, base_list):
    seq_array[np.isin(seq_array, [b"A", b"C", b"G", b"T"], invert=True)] = b"N"
    seq_str_array = [seq.decode('utf-8') for seq in seq_array]
    df = pl.DataFrame({"sequence": seq_str_array})
    df_onehot = df.select(pl.col("sequence").cast(pl.Categorical)).to_dummies()
    for base in base_list:
        base_col = f"sequence_{base.decode()}"
        if base_col not in df_onehot.columns:
            df_onehot = df_onehot.with_column(pl.lit(0).alias(base_col))
    df_onehot = df_onehot.select(sorted(df_onehot.columns))
    return df_onehot.to_numpy()

def parse_encode_list(encode_spec):
    if not encode_spec:
        encode_spec = [b"A", b"C", b"G", b"T", b"N"]
    elif isinstance(encode_spec, (list, tuple)):
        encode_spec = [base.encode() if isinstance(base, str) else base for base in encode_spec]
    elif isinstance(encode_spec, str):
        encode_spec = [base.encode() for base in list(encode_spec)]
    else:
        raise TypeError("Please input string or list of strings!")
    return encode_spec

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
    encode_spec = parse_encode_list(encode_spec)
    return array_to_onehot(seq_data, encode_spec)

def load_chromosome(fasta_file, chrom, encode_spec, output_dir):
    logger.info(f"Encoding chromosome {chrom} from FASTA file {fasta_file}")
    try:
        with FastaFile(fasta_file) as fasta:
            fasta_seq = fasta.fetch(chrom)
            onehot_seq = encode_sequence(fasta_seq, encode_spec)
        tmp_h5_file = os.path.join(output_dir, f'{chrom}.h5')
        with h5py.File(tmp_h5_file, 'w') as f:
            f.create_dataset('sequence', data=onehot_seq, compression=32001, compression_opts=(0, 2, 0, 0, 5, 1, 2))
        logger.info(f"Finished encoding and saving chromosome {chrom} to {tmp_h5_file}")
        return chrom, tmp_h5_file
    except Exception as e:
        logger.error(f"Error encoding chromosome {chrom}: {e}")
        raise

class ReferenceGenome:
    def __init__(self, fasta_file=None, encode_spec=None, hdf5_file=None, output_dir=None):
        self.encode_spec = parse_encode_list(encode_spec)
        self.output_dir = output_dir
        if fasta_file:
            self.fasta_file = fasta_file
            self.genome_df = self.load_genome_parallel()
        elif hdf5_file:
            self.hdf5_file = hdf5_file
            self.genome_df = self.load_from_hdf5(hdf5_file)
        else:
            raise ValueError("Either fasta_file or hdf5_file must be provided.")

    def load_genome_parallel(self):
        chrom_list = [f'chr{i}' for i in range(1, 23)]
        logger.info("Starting parallel encoding of genome")
        with Pool(cpu_count()) as pool:
            results = pool.starmap(load_chromosome, [(self.fasta_file, chrom, self.encode_spec, self.output_dir) for chrom in chrom_list])
        genome_data = [(chrom, h5_file) for chrom, h5_file in results]
        genome_df = pl.DataFrame(genome_data, schema=[("chrom", pl.Utf8), ("file_path", pl.Utf8)])
        logger.info("Finished parallel encoding of genome")
        return genome_df

    def get_sequence(self, chrom, start, end):
        file_path = self.genome_df.filter(pl.col("chrom") == chrom).select(pl.col("file_path"))[0, 0]
        with h5py.File(file_path, 'r') as f:
            sequence = f['sequence'][start:end]
        return np.array(sequence, dtype=np.int8)

    def save_to_hdf5(self, hdf5_file):
        logger.info(f"Saving entire reference genome to {hdf5_file}")
        try:
            with h5py.File(hdf5_file, 'w') as f:
                for row in self.genome_df.iterrows():
                    chrom, file_path = row
                    with h5py.File(file_path, 'r') as tmp_f:
                        seq_data = tmp_f['sequence'][:]
                        seq_group = f.create_group(chrom)
                        seq_group.create_dataset('sequence', data=seq_data, chunks=True, compression=32001, compression_opts=(2, 2, 0, 0, 5, 1, 2))
            logger.info(f"Successfully saved reference genome to {hdf5_file}")
        except Exception as e:
            logger.error(f"Error saving reference genome to {hdf5_file}: {e}")
            raise

    def load_from_hdf5(self, hdf5_file):
        logger.info(f"Loading reference genome from {hdf5_file}")
        try:
            genome_data = []
            with h5py.File(hdf5_file, 'r') as f:
                for chrom in f.keys():
                    sequence = f[chrom]['sequence'][:]
                    genome_data.append((chrom, sequence.tolist()))
            genome_df = pl.DataFrame(genome_data, schema=[("chrom", pl.Utf8), ("sequence", pl.List(pl.Int8))])
            logger.info(f"Successfully loaded reference genome from {hdf5_file}")
            return genome_df
        except Exception as e:
            logger.error(f"Error loading reference genome from {hdf5_file}: {e}")
            raise

def merge_h5_files(tmp_dir, final_h5_file):
    logger.info(f"Merging HDF5 files from {tmp_dir} to {final_h5_file}")
    try:
        with h5py.File(final_h5_file, 'a') as final_file:
            for tmp_file in os.listdir(tmp_dir):
                if tmp_file.endswith(".h5"):
                    tmp_file_path = os.path.join(tmp_dir, tmp_file)
                    logger.info(f"Processing temporary file: {tmp_file_path}")
                    
                    try:
                        with h5py.File(tmp_file_path, 'r') as tmp:
                            for chrom in tmp.keys():
                                logger.info(f"Processing chromosome: {chrom}")
                                
                                # Check if the chromosome is a dataset
                                if isinstance(tmp[chrom], h5py.Dataset):
                                    chrom_group = final_file.require_group(chrom)
                                    if 'sequence' in chrom_group:
                                        del chrom_group['sequence']
                                    
                                    try:
                                        tmp.copy(chrom, chrom_group, name='sequence')
                                    except Exception as e:
                                        logger.error(f"Error copying dataset for chromosome {chrom}: {e}")
                                        raise
                                else:
                                    logger.warning(f"Skipping chromosome {chrom} in {tmp_file_path} as it is not a dataset")
                    
                    except Exception as e:
                        logger.error(f"Error opening temporary file {tmp_file_path}: {e}")
                        raise
        
        logger.info("Finished merging HDF5 files")
    
    except Exception as e:
        logger.error(f"An error occurred while merging HDF5 files: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to encode reference genome and merge chromosome HDF5 files')
    parser.add_argument("--fasta", required=True, type=str, help="Path to reference genome FASTA file")
    parser.add_argument("--outdir", required=True, type=str, help="Path to results save folder")

    args = parser.parse_args()

    output_dir = os.path.join(args.outdir, 'tmp_chrom_files')
    os.makedirs(output_dir, exist_ok=True)

    ref_hdf5_file = os.path.join(args.outdir, 'reference_genome.h5')
    
    try:
        ref_genome = ReferenceGenome(fasta_file=args.fasta, output_dir=output_dir)
        merge_h5_files(output_dir, ref_hdf5_file)
    except Exception as e:
        logger.error(f"An error occurred during processing: {e}")
    finally:
        shutil.rmtree(output_dir)
        logger.info("Cleaned up temporary files")

    logger.info(f"Reference genome HDF5 file created at {ref_hdf5_file}")