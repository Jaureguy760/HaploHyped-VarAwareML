import os
import argparse
import h5py
import polars as pl
import numpy as np
from multiprocessing import Pool, cpu_count
import gzip
import shutil
from utils.common_utils import nucleotide_to_index, bitpack_indices
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,  # Set the logging level
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Log message format
    handlers=[
        logging.FileHandler("haplohyped.log"),  # Log to a file
        logging.StreamHandler()  # Log to the console
    ]
)

logger = logging.getLogger(__name__)

def read_sample_list(sample_list_path):
    """
    Read the sample list file and return a list of sample names.

    Parameters:
    sample_list_path (str): Path to the sample list file.

    Returns:
    list: A list of sample names.
    """
    try:
        with open(sample_list_path, 'r') as f:
            return [line.strip() for line in f]
    except FileNotFoundError as e:
        logger.error(f"Sample list file not found: {e}")
        raise

def store_individuals(data_path, save_path, study_name):
    """
    Store individual sample names into an HDF5 file with compression.

    Parameters:
    data_path (str): Path to the sample list file.
    save_path (str): Path to save the HDF5 file.
    study_name (str): Name of the study.
    """
    logger.info(f"Storing individuals from {data_path} to {save_path}")
    h5_file = os.path.join(save_path, f'{study_name}.h5')
    individuals = read_sample_list(data_path)
    
    with h5py.File(h5_file, 'a') as h5_gen_file:
        if 'individuals' in h5_gen_file:
            del h5_gen_file['individuals']
        h5_gen_file.create_dataset('individuals', data=np.array(individuals, dtype='S'), compression='gzip', compression_opts=5)
    logger.info("Finished storing individuals")

def genotype_VCF2hdf5(data_path, donor_id, chromosome, save_path, study_name):
    """
    Process VCF file and store genotype data into an HDF5 file with compression.

    Parameters:
    data_path (str): Path to the input VCF file.
    donor_id (str): Identifier for the donor.
    chromosome (int): Chromosome number.
    save_path (str): Path to save the HDF5 file.
    study_name (str): Name of the study.
    """
    logger.info(f"Processing VCF file {data_path} for donor {donor_id} and chromosome {chromosome}")
    tmp_h5_file = os.path.join(save_path, f'{study_name}_tmp_{donor_id}_{chromosome}.h5')

    # Read column names from VCF header
    columns = []
    with gzip.open(data_path, 'rt') as file:
        for line in file:
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                break

    # Read VCF data into DataFrame
    df = pl.read_csv(data_path, sep='\t', comment='#', compression='gzip', has_header=False)
    df.columns = columns
    genotype_cols = columns[9:]

    # Filter and transform genotype data
    df_genotypes = df.select(genotype_cols).apply(lambda x: x if isinstance(x, str) and '|' in x else None)
    df_positions = df.select(['#CHROM', 'POS', 'ID', 'REF', 'ALT'])
    
    # Convert nucleotides to indices
    df_positions = df_positions.with_columns([
        pl.col('REF').apply(nucleotide_to_index).alias('REF_IDX'),
        pl.col('ALT').apply(nucleotide_to_index).alias('ALT_IDX')
    ])
    merged_df = pl.concat([df_positions, df_genotypes], how='horizontal')

    # Store genotype data into HDF5
    with h5py.File(tmp_h5_file, 'a') as h5_gen_file:
        group_path = f'donor_{donor_id}/chr_{chromosome}'
        group = h5_gen_file.require_group(group_path)
        if 'genotype' in group:
            del group['genotype']
        packed_genotypes = bitpack_indices(merged_df.select(genotype_cols).to_numpy().astype(np.int8))
        group.create_dataset('genotype', data=packed_genotypes, compression='lzf', chunks=True)
    logger.info(f"Finished processing VCF file for donor {donor_id} and chromosome {chromosome}")

def process_chromosome(chromosome, vcf_dir, out_dir, study_name, donor_ids):
    """
    Process a single chromosome's VCF file and store the data into HDF5 files.

    Parameters:
    chromosome (int): Chromosome number to process.
    vcf_dir (str): Directory containing VCF files.
    out_dir (str): Directory to save the output HDF5 files.
    study_name (str): Name of the study.
    donor_ids (list): List of donor IDs to process.
    """
    logger.info(f"Processing chromosome {chromosome} for donors {donor_ids}")
    vcf_file = os.path.join(vcf_dir, f'chr{chromosome}.filtered.vcf.gz')
    for donor_id in donor_ids:
        genotype_VCF2hdf5(vcf_file, donor_id, chromosome, out_dir, study_name)

def merge_h5_files(tmp_dir, final_h5_file):
    """
    Merge temporary HDF5 files into a single final HDF5 file.

    Parameters:
    tmp_dir (str): Directory containing temporary HDF5 files.
    final_h5_file (str): Path to the final HDF5 file.
    """
    logger.info(f"Merging HDF5 files from {tmp_dir} to {final_h5_file}")
    with h5py.File(final_h5_file, 'a') as final_file:
        for tmp_file in os.listdir(tmp_dir):
            if tmp_file.endswith(".h5"):
                tmp_file_path = os.path.join(tmp_dir, tmp_file)
                with h5py.File(tmp_file_path, 'r') as tmp:
                    for donor in tmp.keys():
                        if donor not in final_file:
                            final_file.create_group(donor)
                        for chrom in tmp[donor].keys():
                            group_path = f"{donor}/{chrom}"
                            if group_path in final_file:
                                del final_file[group_path]
                            tmp.copy(f"{donor}/{chrom}", final_file)
    logger.info("Finished merging HDF5 files")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to convert VCF data')
    parser.add_argument("--cohort_name", required=True, type=str, help="Cohort specific name")
    parser.add_argument("--vcf", required=True, type=str, help="Path to VCF files directory")
    parser.add_argument("--outdir", required=True, type=str, help="Path to results save folder")
    parser.add_argument("--flag", required=True, type=str, choices=['individuals', 'chunk'], help="Type of data to process")
    parser.add_argument("--sample_list", required=False, type=str, help="Path to sample list file", default='sample_list.txt')

    args = parser.parse_args()

    tmp_dir = os.path.join(args.outdir, 'tmp_files')
    os.makedirs(tmp_dir, exist_ok=True)
    
    if args.flag == 'individuals':
        store_individuals(args.sample_list, args.outdir, args.cohort_name)
    elif args.flag == 'chunk':
        donor_ids = read_sample_list(args.sample_list)
        chromosomes = range(1, 23)
        num_processes = cpu_count() - 4
        with Pool(processes=num_processes) as pool:
            pool.starmap(process_chromosome, [(chr, args.vcf, tmp_dir, args.cohort_name, donor_ids) for chr in chromosomes])

        final_h5_file = os.path.join(args.outdir, f'{args.cohort_name}.h5')
        merge_h5_files(tmp_dir, final_h5_file)
        shutil.rmtree(tmp_dir)
