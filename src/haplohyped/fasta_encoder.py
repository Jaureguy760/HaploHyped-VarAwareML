import os
import h5py
import numpy as np
from concurrent.futures import ThreadPoolExecutor
import logging
import polars as pl
from pysam import FastaFile
import shutil
import click
import hdf5plugin
import b2h5py.auto

# Set up logging configuration for both file and console outputs
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("reference_genome.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class ReferenceGenome:
    def __init__(self, fasta_file=None, encode_spec=None, hdf5_file=None, output_dir=None):
        self.encode_spec = self.parse_encode_list(encode_spec)
        self.output_dir = output_dir
        self.fasta_file = fasta_file
        self.hdf5_file = hdf5_file
        self.genome_df = None

    @staticmethod
    def parse_encode_list(encode_spec):
        """
        Parse the encoding specification into a list of bytes.
        """
        if not encode_spec:
            encode_spec = [b"A", b"C", b"G", b"T", b"N"]
        elif isinstance(encode_spec, (list, tuple)):
            encode_spec = [base.encode() if isinstance(base, str) else base for base in encode_spec]
        elif isinstance(encode_spec, str):
            encode_spec = [base.encode() for base in list(encode_spec)]
        else:
            raise TypeError("Please input string or list of strings!")
        return encode_spec

    @staticmethod
    def array_to_onehot(seq_array, base_list):
        """
        Convert a sequence array to a one-hot encoded numpy array.
        """
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

    def encode_sequence(self, seq_data, ignore_case=True):
        """
        Encode the sequence data into a one-hot encoded format.
        """
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
        return self.array_to_onehot(seq_data, self.encode_spec)

    def load_chromosome(self, chrom):
        """
        Load and encode a specific chromosome from the FASTA file.
        """
        logger.info(f"Encoding chromosome {chrom} from FASTA file {self.fasta_file}")
        try:
            with FastaFile(self.fasta_file) as fasta:
                fasta_seq = fasta.fetch(chrom)
                onehot_seq = self.encode_sequence(fasta_seq)
            tmp_h5_file = os.path.join(self.output_dir, f'{chrom}.h5')
            with h5py.File(tmp_h5_file, 'w') as f:
                f.create_dataset('sequence', data=onehot_seq, compression=32001, compression_opts=(0, 2, 0, 0, 5, 1, 2))
            logger.info(f"Finished encoding and saving chromosome {chrom} to {tmp_h5_file}")
            return chrom, tmp_h5_file
        except Exception as e:
            logger.error(f"Error encoding chromosome {chrom}: {e}")
            raise

    def load_genome_parallel(self):
        """
        Load and encode the entire genome in parallel using multiple threads.
        """
        chrom_list = [f'chr{i}' for i in range(1, 23)]
        logger.info("Starting parallel encoding of genome")
        with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
            results = list(executor.map(self.load_chromosome, chrom_list))
        genome_data = [(chrom, h5_file) for chrom, h5_file in results]
        self.genome_df = pl.DataFrame(genome_data, schema=[("chrom", pl.Utf8), ("file_path", pl.Utf8)])
        logger.info("Finished parallel encoding of genome")
        return self.genome_df

    def get_sequence(self, chrom, start, end):
        """
        Retrieve a specific sequence segment from the encoded genome.
        """
        file_path = self.genome_df.filter(pl.col("chrom") == chrom).select(pl.col("file_path"))[0, 0]
        with h5py.File(file_path, 'r') as f:
            sequence = f['sequence'][start:end]
        return np.array(sequence, dtype=np.int8)

class HDF5Handler:
    @staticmethod
    def save_to_hdf5(genome_df, hdf5_file):
        """
        Save the entire reference genome data to a single HDF5 file.
        """
        logger.info(f"Saving entire reference genome to {hdf5_file}")
        try:
            with h5py.File(hdf5_file, 'w') as f:
                for row in genome_df.iterrows():
                    chrom, file_path = row
                    with h5py.File(file_path, 'r') as tmp_f:
                        seq_data = tmp_f['sequence'][:]
                        seq_group = f.create_group(chrom)
                        seq_group.create_dataset('sequence', data=seq_data, chunks=True, compression=32001, compression_opts=(2, 2, 0, 0, 5, 1, 2))
            logger.info(f"Successfully saved reference genome to {hdf5_file}")
        except Exception as e:
            logger.error(f"Error saving reference genome to {hdf5_file}: {e}")
            raise

    @staticmethod
    def load_from_hdf5(hdf5_file):
        """
        Load the reference genome data from a single HDF5 file.
        """
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

    @staticmethod
    def merge_h5_files(tmp_dir, final_h5_file):
        """
        Merge individual HDF5 chromosome files into a single HDF5 file.
        """
        logger.info(f"Merging HDF5 files from {tmp_dir} to {final_h5_file}")
        try:
            with h5py.File(final_h5_file, 'a') as final_file:
                def process_tmp_file(tmp_file):
                    tmp_file_path = os.path.join(tmp_dir, tmp_file)
                    if not tmp_file.endswith(".h5"):
                        return
                    logger.info(f"Processing temporary file: {tmp_file_path}")
                    with h5py.File(tmp_file_path, 'r') as tmp:
                        for chrom in tmp.keys():
                            logger.info(f"Processing chromosome: {chrom}")
                            if isinstance(tmp[chrom], h5py.Dataset):
                                chrom_group = final_file.require_group(chrom)
                                if 'sequence' in chrom_group:
                                    del chrom_group['sequence']
                                tmp.copy(chrom, chrom_group, name='sequence')

                with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
                    executor.map(process_tmp_file, os.listdir(tmp_dir))

            logger.info("Finished merging HDF5 files")
        except Exception as e:
            logger.error(f"An error occurred while merging HDF5 files: {e}")
            raise

@click.command()
@click.option("--fasta", required=True, type=click.Path(exists=True), help="Path to reference genome FASTA file")
@click.option("--outdir", required=True, type=click.Path(), help="Path to results save folder")
@click.option("--cores", default=os.cpu_count(), type=int, help="Number of CPU cores to use")
def main(fasta, outdir, cores):
    """
    Main function to orchestrate the encoding and saving of the reference genome.
    """
    output_dir = os.path.join(outdir, 'tmp_chrom_files')
    os.makedirs(output_dir, exist_ok=True)

    ref_hdf5_file = os.path.join(outdir, 'reference_genome.h5')

    try:
        ref_genome = ReferenceGenome(fasta_file=fasta, output_dir=output_dir)
        ref_genome.load_genome_parallel()
        HDF5Handler.merge_h5_files(output_dir, ref_hdf5_file)
    except Exception as e:
        logger.error(f"An error occurred during processing: {e}")
    finally:
        shutil.rmtree(output_dir)
        logger.info("Cleaned up temporary files")

    logger.info(f"Reference genome HDF5 file created at {ref_hdf5_file}")

if __name__ == "__main__":
    main()
