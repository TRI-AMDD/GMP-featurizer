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
        int calculate_solid_gmpordernorm_noderiv_elemental_sigma_cutoff(double**, double**, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**);

        int calculate_solid_gmpordernorm_noderiv_sigma_cutoff(double**, double**, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);
        
        int calculate_solid_gmpordernorm_noderiv_original(double**, double**, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);
        
        int calculate_gmpordernorm_noderiv(double**, double**, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);

        int calculate_solid_gmpordernorm(double**, double**, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);
        
        int calculate_gmpordernorm(double**, double**, double**, double**, double**, int*,
                                    int*, int, int,
                                    int**, double**, int, double**, int*, int*,
                                    double**, double**);
    """
)
ffibuilder.set_source(
    "GMPFeaturizer.GMPOrderNorm._libgmpordernorm",
    '#include "calculate_gmpordernorm.h"',
    sources=[
        "GMPFeaturizer/GMPOrderNorm/calculate_gmpordernorm.cpp",
        # "GMP-featurizer/GMPOrderNorm/gmpordernorm.cpp",
        "GMPFeaturizer/GMPOrderNorm/helper.cpp",
        "GMPFeaturizer/GMPOrderNorm/surface_harmonics.cpp",
        "GMPFeaturizer/GMPOrderNorm/solid_harmonics.cpp",
    ],
    source_extension=".cpp",
    include_dirs=["GMPFeaturizer/GMPOrderNorm/"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
