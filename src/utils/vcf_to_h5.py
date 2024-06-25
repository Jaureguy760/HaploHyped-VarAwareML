# VCFtoHDF5Converter/vcf_to_hdf5.py

import os
import parse_vcf
import polars as pl
import h5py
import numpy as np
import logging
import time
import shutil
import hdf5plugin
import b2h5py.auto
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import List
import click

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("haplohyped.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class VCFtoHDF5Converter:
    """
    Class to convert VCF files to HDF5 format.
    """

    def __init__(self, cohort_name: str, vcf_dir: str, out_dir: str, sample_list_path: str, cores: int, cxx_threads: int):
        """
        Initialize the converter with necessary parameters.

        Args:
            cohort_name (str): Cohort specific name.
            vcf_dir (str): Path to VCF files directory.
            out_dir (str): Path to results save folder.
            sample_list_path (str): Path to sample list file.
            cores (int): Number of CPU cores to use.
            cxx_threads (int): Number of threads to use in the C++ code.
        """
        self.cohort_name = cohort_name
        self.vcf_dir = vcf_dir
        self.out_dir = out_dir
        self.sample_list_path = sample_list_path
        self.cores = cores
        self.cxx_threads = cxx_threads
        self.donor_ids = self.read_sample_list(sample_list_path)
        self.chromosomes = range(1, 23)
        self.tmp_dir = os.path.join(out_dir, 'tmp_files')
        os.makedirs(self.tmp_dir, exist_ok=True)

    def read_sample_list(self, sample_list_path: str) -> List[str]:
        """
        Read and return a list of sample IDs from a file.

        Args:
            sample_list_path (str): Path to the sample list file.

        Returns:
            List[str]: List of sample IDs.

        Raises:
            FileNotFoundError: If the sample list file is not found.
            Exception: For any other error during file reading.
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

    def genotype_VCF2hdf5(self, data_path: str, donor_id: str, chromosome: int) -> None:
        """
        Convert VCF data to HDF5 format for a specific donor and chromosome.

        Args:
            data_path (str): Path to the VCF file.
            donor_id (str): ID of the donor.
            chromosome (int): Chromosome number.

        Raises:
            Exception: For any error during VCF processing or HDF5 file creation.
        """
        # Set the number of threads for C++ code
        os.environ["OMP_NUM_THREADS"] = str(self.cxx_threads)
        os.environ["MKL_NUM_THREADS"] = str(self.cxx_threads)
        
        logger.info(f"Processing VCF file {data_path} for chromosome {chromosome} and donor {donor_id}")
        tmp_h5_file = os.path.join(self.tmp_dir, f'{self.cohort_name}_tmp_donor_{donor_id}_chr_{chromosome}.h5')

        try:
            chrom_str = f"chr{chromosome}"
            
            if donor_id:
                # Load VCF data using the parse_vcf module
                results = parse_vcf.load_vcf(data_path, donor_id, chrom_str)
                columns = ["chrom", "start", "stop", "ref", "alt", "phase1", "phase2"]
                snp_df = pl.DataFrame({
                    "chrom": [r[0] for r in results],
                    "start": [r[1] for r in results],
                    "stop": [r[2] for r in results],
                    "ref": [r[3] for r in results],
                    "alt": [r[4] for r in results],
                    "phase1": [r[5] for r in results],
                    "phase2": [r[6] for r in results]
                })
                snp_df = snp_df.with_columns([
                    pl.col("start").cast(pl.UInt32),
                    pl.col("stop").cast(pl.UInt32),
                    pl.col("phase1").cast(pl.Int8),
                    pl.col("phase2").cast(pl.Int8)
                ])

                dtype = np.dtype([
                    ('chrom', 'S5'),
                    ('start', np.uint32),
                    ('stop', np.uint32),
                    ('ref', 'S10'),
                    ('alt', 'S10'),
                    ('phase1', np.int8),
                    ('phase2', np.int8)
                ])

                snp_struct = np.array([tuple(row) for row in snp_df.iter_rows()], dtype=dtype)

                # Write data to HDF5 file
                with h5py.File(tmp_h5_file, 'w') as h5_gen_file:
                    group_path = f'donor_{donor_id}/chr_{chromosome}'
                    group = h5_gen_file.create_group(group_path)
                    group.create_dataset('snp_data', data=snp_struct, compression=32001,
                                        compression_opts=(2, 2, 0, 0, 5, 1, 2), chunks=True)

                logger.info(f"Finished processing VCF file for donor {donor_id} and chromosome {chromosome}")
        except Exception as e:
            logger.error(f"An error occurred while processing VCF file: {e}")
            raise

    def process_donor(self, donor_id: str) -> None:
        """
        Process all chromosomes for a single donor.

        Args:
            donor_id (str): ID of the donor.
        """
        logger.info(f"Processing donor {donor_id}")
        for chromosome in self.chromosomes:
            vcf_file = os.path.join(self.vcf_dir, f'chr{chromosome}.filtered.vcf.gz')
            self.genotype_VCF2hdf5(vcf_file, donor_id, chromosome)

    def merge_h5_files(self) -> None:
        """
        Merge multiple HDF5 files into a single file.

        Raises:
            Exception: For any error during the merging process.
        """
        final_h5_file = os.path.join(self.out_dir, f'{self.cohort_name}.h5')
        logger.info(f"Merging HDF5 files from {self.tmp_dir} to {final_h5_file}")
        try:
            with h5py.File(final_h5_file, 'w') as final_file:
                for tmp_file in os.listdir(self.tmp_dir):
                    if tmp_file.endswith(".h5"):
                        tmp_file_path = os.path.join(self.tmp_dir, tmp_file)
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

    def run(self):
        """
        Run the VCF to HDF5 conversion process.
        """
        start_time = time.time()

        try:
            process_donor_partial = partial(self.process_donor)

            with ProcessPoolExecutor(max_workers=self.cores) as executor:
                executor.map(process_donor_partial, self.donor_ids)

            merge_start_time = time.time()
            self.merge_h5_files()
            merge_end_time = time.time()

            end_time = time.time()
            total_time = end_time - start_time

            logger.info(f"Time taken to merge HDF5 files: {merge_end_time - merge_start_time:.2f} seconds")
            logger.info(f"Total time taken: {total_time:.2f} seconds")

        except Exception as e:
            logger.error(f"An error occurred: {e}")
        finally:
            shutil.rmtree(self.tmp_dir)

@click.command()
@click.option('--cohort_name', required=True, type=str, help="Cohort specific name")
@click.option('--vcf', required=True, type=str, help="Path to VCF files directory")
@click.option('--outdir', required=True, type=str, help="Path to results save folder")
@click.option('--sample_list', required=True, type=str, help="Path to sample list file")
@click.option('--cores', default=os.cpu_count(), type=int, help="Number of CPU cores to use")
@click.option('--cxx_threads', default=4, type=int, help="Number of threads to use in the C++ code")
def main(cohort_name, vcf, outdir, sample_list, cores, cxx_threads):
    """
    Main function to parse arguments and run the converter.
    """
    converter = VCFtoHDF5Converter(
        cohort_name=cohort_name,
        vcf_dir=vcf,
        out_dir=outdir,
        sample_list_path=sample_list,
        cores=cores,
        cxx_threads=cxx_threads
    )

    converter.run()

if __name__ == "__main__":
    main()
