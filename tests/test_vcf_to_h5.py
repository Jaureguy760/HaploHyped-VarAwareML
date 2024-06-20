# scripts/test_benchmark.py

import os, sys
import time
import h5py
import numpy as np
import logging
from subprocess import call
import shutil
# sys.path.append('/iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML/src/utils')

from .common_utils import nucleotide_to_index, bitpack_indices
# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("haplohyped_test.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_command(command):
    """
    Run a shell command and log the output.

    Parameters:
    command (str): The command to run.
    """
    logger.info(f"Running command: {command}")
    start_time = time.time()
    result = call(command, shell=True)
    end_time = time.time()
    if result != 0:
        logger.error(f"Command failed: {command}")
    else:
        logger.info(f"Command succeeded: {command} (Time: {end_time - start_time:.2f} seconds)")

def verify_h5_file(h5_file):
    """
    Verify the contents of the HDF5 file.

    Parameters:
    h5_file (str): Path to the HDF5 file.
    """
    try:
        with h5py.File(h5_file, 'r') as f:
            logger.info(f"Verifying HDF5 file: {h5_file}")
            for donor in f.keys():
                for chrom in f[donor].keys():
                    data = f[f"{donor}/{chrom}/genotype"]
                    if data.shape[0] == 0:
                        logger.warning(f"No data found for {donor}/{chrom}")
                    else:
                        logger.info(f"Data verified for {donor}/{chrom}")
    except Exception as e:
        logger.error(f"An error occurred while verifying HDF5 file: {e}")

def benchmark_h5_file(h5_file):
    """
    Benchmark the I/O performance of the HDF5 file.

    Parameters:
    h5_file (str): Path to the HDF5 file.
    """
    try:
        with h5py.File(h5_file, 'r') as f:
            logger.info(f"Benchmarking HDF5 file: {h5_file}")
            start_time = time.time()
            for donor in f.keys():
                for chrom in f[donor].keys():
                    data = f[f"{donor}/{chrom}/genotype"][...]
                    _ = np.mean(data)  # Perform a dummy operation
            end_time = time.time()
            logger.info(f"Benchmark completed: {h5_file} (Time: {end_time - start_time:.2f} seconds)")
    except Exception as e:
        logger.error(f"An error occurred during benchmarking: {e}")

def main():
    cohort_name = "ipscs"
    vcf_dir = "/iblm/netapp/data4/Frazer_collab/ipscs/datasets/raw/genotypes/michigan_impute/results/final_vcf/merged_sorted.vcf.gz"
    out_dir = "/iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML/tests/out"
    sample_list = "/iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML/tests/data/ipscs_samples_test.txt"
    
    # Prepare temporary directory
    tmp_dir = os.path.join(out_dir, 'tmp_files')
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Store individuals
    # command = f"python process_vcf.py --cohort_name {cohort_name} --vcf {vcf_dir} --outdir {out_dir} --flag individuals --sample_list {sample_list}"
    # run_command(command)
    
    # Process VCF to HDF5 (Chunk-wise)
    command = f"python /iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML/src/utils/vcf_to_h5.py --cohort_name {cohort_name} --vcf {vcf_dir} --outdir {out_dir} --flag chunk --sample_list {sample_list}"
    run_command(command)
    
    # Verify the output HDF5 file
    final_h5_file = os.path.join(out_dir, f'{cohort_name}.h5')
    verify_h5_file(final_h5_file)
    
    # Benchmark the output HDF5 file
    benchmark_h5_file(final_h5_file)
    
    # Clean up temporary directory
    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    main()
