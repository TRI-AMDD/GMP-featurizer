#include <math.h>
//#include "mpi.h"
// #include "gmpordernorm.h"
// #include "helper.h"
#include "surface_harmonics.h"
#include "solid_harmonics.h"
#include "solid_harmonics_optimization.h"
#include "solid_harmonics_optimization2.h"
#include "solid_harmonics_optimization3.h"

        
extern "C" int calculate_surface_gmpordernorm_fp_deriv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

extern "C" int calculate_surface_gmpordernorm_noderiv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);



extern "C" int calculate_solid_gmpordernorm_noderiv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_occ_deriv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_fp_deriv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_fp_occ_deriv_ref(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**, double**);



extern "C" int calculate_solid_gmpordernorm_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_fp_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_fp_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**, double**);



extern "C" int calculate_solid_gmpordernorm_sigma_cutoff_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_sigma_cutoff_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_sigma_cutoff_fp_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_sigma_cutoff_fp_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int,
                                        int**, double**, int, double**, int*, int*,
                                        double**, double**, double**);



extern "C" int calculate_solid_gmpordernorm_elemental_sigma_cutoff_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_cutoff_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int,
                                        int**, double**, int, double**, int*, double**, int*,
                                        double**, double**, double**);



extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double*, double*, double*, double*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt2(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double*, double*, double*, double*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt2_2(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double*, double*, double*, double*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt3(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double*, double*, double*, double*,
                                        double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt3_2(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double*, double*, double*, double*,
                                        double**);


extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_occ_deriv_opt2_2(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double*, double*, double*, double*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double**, double**);

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_occ_deriv(double**, double**, double*, double**, double**, double**, int*,
                                        int*, int, int, int, int,
                                        int**, double**, int, double**, int*, double**, double**, int*,
                                        double**, double**, double**);





bool check_implementation(int nmcsh,int** params_i);
void calculate_bin_ranges(double** cell, double** scale, int natoms, double cutoff, int& max_atoms_bin, 
                          int& neigh_check_bins,int **bin_i, int* bin_range, int* nbins);
void find_neighbors(double** cell, double** cart, double* occupancies, double** ref_cart, 
                    double** scale, double** ref_scale, int* pbc_bools,int* atom_i, int natoms,
                    int calc_ii, double cutoff_sqr, int* bin_range, int* nbins, int** bin_i,
                    int& nneigh, double* nei_list_d, int* nei_list_i, double* nei_list_occupancy);

const int NUM_IMPLEMENTED_TYPE = 73;
const int IMPLEMENTED_MCSH_TYPE[][2] = {
    {-1,0},
    {0, 0},
    {1, 0},
    {2, 0},
    {3, 0},
    {4, 0},
    {5, 0},
    {6, 0},
    {7, 0},
    {8, 0},
    {9, 0},
    {-1,1},
    {0, 1},
    {1, 1},
    {2, 1},
    {3, 1},
    {4, 1},
    {5, 1},
    {6, 1},
    {7, 1},
    {8, 1},
    {9, 1},
    {0, 1},
    {1, 1},
    {2, 1},
    {2, 2},
    {3, 1},
    {3, 2},
    {3, 3},
    {4, 1},
    {4, 2},
    {4, 3},
    {4, 4},
    {5, 1},
    {5, 2},
    {5, 3},
    {5, 4},
    {5, 5},
    {6, 1},
    {6, 2},
    {6, 3},
    {6, 4},
    {6, 5},
    {6, 6},
    {6, 7},
    {7, 1},
    {7, 2},
    {7, 3},
    {7, 4},
    {7, 5},
    {7, 6},
    {7, 7},
    {7, 8},
    {8, 1},
    {8, 2},
    {8, 3},
    {8, 4},
    {8, 5},
    {8, 6},
    {8, 7},
    {8, 8},
    {8, 9},
    {8, 10},
    {9, 1},
    {9, 2},
    {9, 3},
    {9, 4},
    {9, 5},
    {9, 6},
    {9, 7},
    {9, 8},
    {9, 9},
    {9, 10},
    {9, 11},
    {9, 12}
}; // Change this when you implement new type!


