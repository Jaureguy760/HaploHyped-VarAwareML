import gzip
import pandas as pd

def check_vcf_columns(vcf_file):
    # Read the VCF file header
    with gzip.open(vcf_file, 'rt') as file:
        for line in file:
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                print("Columns in VCF file:")
                print(columns)
                break

    # Read the first non-header line to verify the data
    df = pd.read_csv(vcf_file, sep='\t', comment='#', compression='gzip', header=None, nrows=1)
    print("First data row in VCF file:")
    print(df)

vcf_file = "/iblm/netapp/data4/Frazer_collab/ipscs/datasets/raw/genotypes/michigan_impute/results/final_vcf/chrom_split/chr10.filtered.vcf.gz"
check_vcf_columns(vcf_file)
