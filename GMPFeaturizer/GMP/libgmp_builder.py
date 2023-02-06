import cffi

ffibuilder = cffi.FFI()
ffibuilder.cdef(
    """
        int calculate_surface_gmpordernorm_fp_deriv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

        int calculate_surface_gmpordernorm_noderiv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);



        int calculate_solid_gmpordernorm_noderiv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);

        int calculate_solid_gmpordernorm_occ_deriv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_fp_deriv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_fp_occ_deriv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**, double**);



        int calculate_solid_gmpordernorm_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);

        int calculate_solid_gmpordernorm_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_fp_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_fp_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**, double**);



        int calculate_solid_gmpordernorm_sigma_cutoff_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);

        int calculate_solid_gmpordernorm_sigma_cutoff_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_sigma_cutoff_fp_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_sigma_cutoff_fp_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**, double**);



        int calculate_solid_gmpordernorm_elemental_sigma_cutoff_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**);

        int calculate_solid_gmpordernorm_elemental_sigma_cutoff_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**, double**, double**);



        int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double**);

        int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
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
        # "GMP-featurizer/GMPOrderNorm/gmpordernorm.cpp",
        "GMPFeaturizer/GMP/helper.cpp",
        "GMPFeaturizer/GMP/surface_harmonics.cpp",
        "GMPFeaturizer/GMP/solid_harmonics.cpp",
    ],
    source_extension=".cpp",
    include_dirs=["GMPFeaturizer/GMP/"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
