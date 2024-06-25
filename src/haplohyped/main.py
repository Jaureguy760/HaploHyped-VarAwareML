# src/haplohyped/main.py

import click
from haplohyped import vcf_to_hdf5, fasta_encoder

@click.group()
def main():
    """HaploHyped Command Line Interface."""
    pass

main.add_command(vcf_to_hdf5.main, name='vcf_to_hdf5')
main.add_command(fasta_encoder.main, name='fasta_encoder')

if __name__ == "__main__":
    main()
