from setuptools import setup, find_packages

setup(
    name='HaploHyped-VarAwareML',
    version='0.1.0',
    author='Jeff Jaureguy',
    author_email='jjaureguy@ucsd.edu',
    description='High-performance VCF to HDF5 pipeline for genomic machine learning',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Jaureguy760/HaploHyped-VarAwareML',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        'click',
        'numpy',
        'h5py',
        'pysam',
        'polars',
        'hdf5plugin',
        'b2h5py',
    ],
    entry_points={
        'console_scripts': [
            'vcf_to_h5=haplohyped.vcf_to_h5:main',
            'fasta_encoder=haplohyped.fasta_encoder:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.8',
)
