# Copyright Toyota Research Institute 2023
"""
Module for compiling the c++ code
"""

import cffi

mkl_include_dir = "/opt/intel/oneapi/mkl/latest/include/"
mkl_lib_dir = "/opt/intel/oneapi/mkl/latest/lib/intel64/"

ffibuilder = cffi.FFI()
ffibuilder.cdef(
    """

    int calculate_solid_gmp_elemental_sigma_gaussian_cutoff_noderiv_opt2(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double*, double*, double*, double*,
        double**);
    """
)
ffibuilder.set_source(
    "GMPFeaturizer.GMP_individual._libgmpind",
    '#include "calculate_gmp.h"',
    sources=[
        "GMPFeaturizer/GMP_individual/calculate_gmp.cpp",
        "GMPFeaturizer/GMP_individual/helper.cpp",
        "GMPFeaturizer/GMP_individual/solid_harmonics_optimization2.cpp",
    ],
    # source_extension=".cpp",
    # include_dirs=["GMPFeaturizer/GMP/"],
    # extra_compile_args=["-g", "-O2"],

    source_extension=".cpp",
    include_dirs=["GMPFeaturizer/GMP_individual/", mkl_include_dir],
    extra_compile_args=["-g", "-O2", "-m64", "-I" + mkl_include_dir],
    libraries=["mkl_rt"],  # Use the runtime MKL library, simplifies linking
    library_dirs=[mkl_lib_dir],
    extra_link_args=["-Wl,-rpath," + mkl_lib_dir],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
