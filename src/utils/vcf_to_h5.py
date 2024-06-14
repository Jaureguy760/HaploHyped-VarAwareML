import os
import argparse
import h5py
import pandas as pd
from multiprocessing import Pool, cpu_count
import gzip
import shutil

def read_sample_list(sample_list_path):
    with open(sample_list_path, 'r') as f:
        return [line.strip() for line in f]

def store_individuals(data_path, save_path, study_name):
    h5_file = os.path.join(save_path, f'{study_name}.h5')
    individuals = []
    with open(data_path, 'r') as f:
        individuals = [line.strip() for line in f]
    with h5py.File(h5_file, 'a') as h5_gen_file:
        if 'individuals' in h5_gen_file:
            del h5_gen_file['individuals']
        h5_gen_file.create_dataset('individuals', data=np.array(individuals, dtype='S'), compression='gzip', compression_opts=5)

def genotype_VCF2hdf5(data_path, donor_id, chromosome, save_path, study_name):
    tmp_h5_file = os.path.join(save_path, f'{study_name}_tmp_{donor_id}_{chromosome}.h5')

    columns = []
    with gzip.open(data_path, 'rt') as file:
        for line in file:
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                break

    df = pd.read_csv(data_path, sep='\t', comment='#', compression='gzip', header=None)
    df.columns = columns
    genotype_cols = columns[9:]

    df_genotypes = df[genotype_cols].applymap(lambda x: x if isinstance(x, str) and '|' in x else None)
    df_positions = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT']]
    merged_df = pd.concat([df_positions, df_genotypes], axis=1)

    with h5py.File(tmp_h5_file, 'a') as h5_gen_file:
        group_path = f'donor_{donor_id}/chr_{chromosome}'
        group = h5_gen_file.require_group(group_path)
        if 'genotype' in group:
            del group['genotype']
        group.create_dataset('genotype', data=merged_df.to_numpy(), compression='gzip', compression_opts=5, chunks=True)

def process_chromosome(chromosome, vcf_dir, out_dir, study_name, donor_ids):
    vcf_file = os.path.join(vcf_dir, f'chr{chromosome}.filtered.vcf.gz')
    for donor_id in donor_ids:
        genotype_VCF2hdf5(vcf_file, donor_id, chromosome, out_dir, study_name)

def merge_h5_files(tmp_dir, final_h5_file):
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
