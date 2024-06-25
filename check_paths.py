import os
import sys

def check_directory(path, name):
    if os.path.isdir(path):
        print(f"{name} directory exists: {path}")
    else:
        print(f"Error: {name} directory does not exist: {path}")
        sys.exit(1)

def check_file(path, name):
    if os.path.isfile(path):
        print(f"{name} file exists: {path}")
    else:
        print(f"Error: {name} file does not exist: {path}")
        sys.exit(1)

# Paths to check
conda_prefix = "/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML"
vcfpp_headers = "/iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML/cpp"

paths_to_check = {
    "VCFPP headers": vcfpp_headers,
    "Conda include": os.path.join(conda_prefix, 'include'),
    "Conda HTSlib include": os.path.join(conda_prefix, 'include', 'htslib'),
    "Conda lib": os.path.join(conda_prefix, 'lib'),
    "libhts.a": os.path.join(conda_prefix, 'lib', 'libhts.a')
}

# Check directories and files
for name, path in paths_to_check.items():
    if "libhts.a" in name:
        check_file(path, name)
    else:
        check_directory(path, name)
