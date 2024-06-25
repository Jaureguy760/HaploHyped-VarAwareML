from setuptools import setup, find_packages

setup(
    name='HaploHyped-VarAwareML',
    version='0.1.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='An end-to-end solution for processing genomic data and performing GPU-accelerated haplotype encoding for machine learning.',
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
            'HaploHyped=haplohyped.main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
)
