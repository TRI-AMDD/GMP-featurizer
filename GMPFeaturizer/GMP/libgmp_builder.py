import cffi

# """int calculate_gmpordernorm(double **, double **, double **, int*,
#                         int *, int, int*, int,
#                         int**, double **, int, double **, int*, int*,
#                         double**, double**);

#         int calculate_gmpordernorm_noderiv(double **, double **, double **, int*,
#                                     int *, int, int*, int,
#                                     int**, double **, int, double **, int*, int*,
#                                     double**);

#     """

ffibuilder = cffi.FFI()
ffibuilder.cdef(
    """
        int calculate_solid_gmpordernorm_noderiv_elemental_sigma_cutoff(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**);

        int calculate_solid_gmpordernorm_noderiv_sigma_cutoff(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);
        
        int calculate_solid_gmpordernorm_noderiv_original(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);
        
        int calculate_gmpordernorm_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);

        int calculate_solid_gmpordernorm_elemental_sigma_cutoff(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**, double**);

        int calculate_solid_gmpordernorm(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);
        
        int calculate_gmpordernorm(double**, double**, double*, double**, double**, double**, int*,
                                    int*, int, int,
                                    int**, double**, int, double**, int*, int*,
                                    double**, double**);
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
