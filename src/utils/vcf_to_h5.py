import os
import argparse
import h5py
import numpy as np
from multiprocessing import Pool, cpu_count
import logging
import polars as pl
from pysam import VariantFile
import shutil
import hdf5plugin  # Ensure this is imported to enable the plugins
import b2h5py.auto  # Automatically enable Blosc2

# Configure logging to output both to a file and to the console
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("haplohyped.log"),
        logging.StreamHandler()
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
    except Exception as e:
        logger.error(f"An error occurred while reading the sample list: {e}")
        raise

def load_vcf(in_vcf, chrom=None, sample=None):
    """
    Load VCF data and return a DataFrame with SNP information.

    Parameters:
    in_vcf (str): Path to the input VCF file.
    chrom (str, optional): Chromosome to filter the data. Defaults to None.
    sample (str, optional): Sample name to filter the data. Defaults to None.

    Returns:
    pl.DataFrame: DataFrame containing SNP data.
    """
    if sample:
        with VariantFile(in_vcf, "r") as vcf:
            vcf.subset_samples([sample])
            vcf_data = vcf.fetch(chrom)
            snp_list = [(record.contig, record.start, record.stop,
                         record.samples[sample]['GT'][0], record.samples[sample]['GT'][1])
                        for record in vcf_data if ((len(record.ref) == 1) and (len(record.alts) == 1)
                                                   and (len(record.alts[0]) == 1))]
            snp_df = pl.DataFrame({
                "chrom": [row[0] for row in snp_list],
                "start": [row[1] for row in snp_list],
                "stop": [row[2] for row in snp_list],
                "phase1": [row[3] for row in snp_list],
                "phase2": [row[4] for row in snp_list],
            }).with_columns([
                pl.col("start").cast(pl.UInt32),
                pl.col("stop").cast(pl.UInt32),
                pl.col("phase1").cast(pl.UInt8),
                pl.col("phase2").cast(pl.UInt8)
            ])
    else:
        logger.error("Data is unphased, please use only phased data for now")
        raise ValueError("Data is unphased, please use only phased data for now")

    return snp_df

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
    logger.info(f"Processing VCF file {data_path} for chromosome {chromosome}")
    tmp_h5_file = os.path.join(save_path, f'{study_name}_tmp_donor_{donor_id}_chr_{chromosome}.h5')

    try:
        snp_df = load_vcf(data_path, chrom=f"chr{chromosome}", sample=donor_id)
        # Ensure the correct data types for columns
        snp_df = snp_df.with_columns([
            pl.col("start").cast(pl.UInt32),
            pl.col("stop").cast(pl.UInt32),
            pl.col("phase1").cast(pl.UInt8),
            pl.col("phase2").cast(pl.UInt8)
        ])
        
        # Convert DataFrame to structured NumPy array with appropriate data types
        snp_struct = np.core.records.fromarrays(
            [snp_df['start'].to_numpy(),
             snp_df['stop'].to_numpy(),
             snp_df['phase1'].to_numpy(),
             snp_df['phase2'].to_numpy()],
            names='start,stop,phase1,phase2'
        )
        
        # Store data into HDF5 with Blosc2 compression and auto-chunking
        with h5py.File(tmp_h5_file, 'a') as h5_gen_file:
            group_path = f'donor_{donor_id}/chr_{chromosome}'
            group = h5_gen_file.require_group(group_path)
            if 'snp_data' in group:
                del group['snp_data']
            group.create_dataset('snp_data', data=snp_struct, compression=32001,
                                 compression_opts=(0, 2, 0, 0, 5, 1, 2))  # Blosc2 + Bitshuffle + LZ4
        logger.info(f"Finished processing VCF file for donor {donor_id} and chromosome {chromosome}")
    except Exception as e:
        logger.error(f"An error occurred while processing VCF file: {e}")
        raise

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
    try:
        with h5py.File(final_h5_file, 'a') as final_file:
            for tmp_file in os.listdir(tmp_dir):
                if tmp_file.endswith(".h5"):
                    tmp_file_path = os.path.join(tmp_dir, tmp_file)
                    with h5py.File(tmp_file_path, 'r') as tmp:
                        for donor in tmp.keys():
                            donor_group = final_file.require_group(donor)
                            for chrom in tmp[donor].keys():
                                chrom_group = donor_group.require_group(chrom)
                                for dset_name in tmp[donor][chrom].keys():
                                    if dset_name in chrom_group:
                                        del chrom_group[dset_name]
                                    tmp.copy(f"{donor}/{chrom}/{dset_name}", chrom_group)
        logger.info("Finished merging HDF5 files")
    except Exception as e:
        logger.error(f"An error occurred while merging HDF5 files: {e}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to convert VCF data')
    parser.add_argument("--cohort_name", required=True, type=str, help="Cohort specific name")
    parser.add_argument("--vcf", required=True, type=str, help="Path to VCF files directory")
    parser.add_argument("--outdir", required=True, type=str, help="Path to results save folder")
    parser.add_argument("--sample_list", required=True, type=str, help="Path to sample list file")

    args = parser.parse_args()

    tmp_dir = os.path.join(args.outdir, 'tmp_files')
    os.makedirs(tmp_dir, exist_ok=True)
    
    donor_ids = read_sample_list(args.sample_list)
    chromosomes = range(4, 5)
    num_processes = cpu_count() - 4
    with Pool(processes=num_processes) as pool:
        pool.starmap(process_chromosome, [(chrom, args.vcf, tmp_dir, args.cohort_name, donor_ids) for chrom in chromosomes])

    final_h5_file = os.path.join(args.outdir, f'{args.cohort_name}.h5')
    merge_h5_files(tmp_dir, final_h5_file)
    shutil.rmtree(tmp_dir)
