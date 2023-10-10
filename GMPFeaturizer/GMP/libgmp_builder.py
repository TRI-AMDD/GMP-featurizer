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
    int calculate_surface_gmpordernorm_fp_deriv_ref(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**);

    int calculate_surface_gmpordernorm_noderiv_ref(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**);



    int calculate_solid_gmpordernorm_noderiv_ref(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**);

    int calculate_solid_gmpordernorm_occ_deriv_ref(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_fp_deriv_ref(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_fp_occ_deriv_ref(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**, double**);



    int calculate_solid_gmpordernorm_noderiv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**);

    int calculate_solid_gmpordernorm_occ_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_fp_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_fp_occ_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**, double**);



    int calculate_solid_gmpordernorm_sigma_cutoff_noderiv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**);

    int calculate_solid_gmpordernorm_sigma_cutoff_occ_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_sigma_cutoff_fp_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_sigma_cutoff_fp_occ_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int,
        int**, double**, int, double**, int*, int*,
        double**, double**, double**);



    int calculate_solid_gmpordernorm_elemental_sigma_cutoff_noderiv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int,
        int**, double**, int, double**, int*, double**, int*,
        double**);

    int calculate_solid_gmpordernorm_elemental_sigma_cutoff_occ_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int,
        int**, double**, int, double**, int*, double**, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int,
        int**, double**, int, double**, int*, double**, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_occ_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int,
        int**, double**, int, double**, int*, double**, int*,
        double**, double**, double**);



    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double*, double*, double*, double*,
        double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt2(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double*, double*, double*, double*,
        double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt2_2(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double*, double*, double*, double*,
        double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt3(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double*, double*, double*, double*,
        double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt3_2(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double*, double*, double*, double*,
        double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_occ_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_occ_deriv_opt2_2(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double*, double*, double*, double*,
        double**, double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double**, double**);

    int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_occ_deriv(
        double**, double**, double*, double**, double**, double**, int*,
        int*, int, int, int, int,
        int**, double**, int, double**, int*, double**, double**, int*,
        double**, double**, double**);
    """
)
ffibuilder.set_source(
    "GMPFeaturizer.GMP._libgmp",
    '#include "calculate_gmpordernorm.h"',
    sources=[
        "GMPFeaturizer/GMP/calculate_gmpordernorm.cpp",
        "GMPFeaturizer/GMP/helper.cpp",
        "GMPFeaturizer/GMP/helper_optimization2.cpp",
        "GMPFeaturizer/GMP/surface_harmonics.cpp",
        "GMPFeaturizer/GMP/solid_harmonics.cpp",
        "GMPFeaturizer/GMP/solid_harmonics_optimization.cpp",
        "GMPFeaturizer/GMP/solid_harmonics_optimization2.cpp",
        "GMPFeaturizer/GMP/solid_harmonics_optimization3.cpp",
    ],
    # source_extension=".cpp",
    # include_dirs=["GMPFeaturizer/GMP/"],
    # extra_compile_args=["-g", "-O2"],

    source_extension=".cpp",
    include_dirs=["GMPFeaturizer/GMP/", mkl_include_dir],
    extra_compile_args=["-g", "-O2", "-m64", "-I" + mkl_include_dir],
    libraries=["mkl_rt"],  # Use the runtime MKL library, simplifies linking
    library_dirs=[mkl_lib_dir],
    extra_link_args=["-Wl,-rpath," + mkl_lib_dir],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
