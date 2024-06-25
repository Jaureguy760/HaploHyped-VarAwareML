cmake_minimum_required(VERSION 3.12)
project(parse_vcf)

set(CMAKE_CXX_STANDARD 14)  # Use C++14

# Set the build type to Debug
set(CMAKE_BUILD_TYPE Debug)

# Include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/include)
include_directories(/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/include/htslib)
include_directories(/iblm/netapp/data4/jjaureguy/vcf_cuda/HaploHyped-VarAwareML/cpp)  # This is the path to vcfpp.h

# Link directories
link_directories(/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib)

# Find Python
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

# Find pybind11
find_package(pybind11 REQUIRED)

# Add library for the Python module
pybind11_add_module(parse_vcf parse_vcf.cpp)

# Link necessary libraries
target_link_libraries(parse_vcf PRIVATE ${Python3_LIBRARIES} 
    /iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib/libhts.a 
    /iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib/liblzma.so 
    /iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib/libdeflate.so 
    /iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib/libbz2.so
    /iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib/libz.so
    -Wl,-rpath,/iblm/netapp/home/jjaureguy/mambaforge/envs/HaploHyped-VarAwareML/lib)
