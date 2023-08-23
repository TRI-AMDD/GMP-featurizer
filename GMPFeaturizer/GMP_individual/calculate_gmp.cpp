#include <math.h>
#include <stdio.h>
#include "calculate_gmp.h"
#include <iostream>
#include <memory>
#include <cassert>


// ##################################################
// ####             Helper function              ####
// ##################################################
bool check_implementation(int nmcsh,int** params_i){
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return false;
    }
    return true;
}

void calculate_bin_ranges(double** cell, double** scale, int natoms, double cutoff, int& max_atoms_bin, 
                          int& neigh_check_bins,int **bin_i, int* bin_range, int* nbins){
    int total_bins;
    double vol, tmp;
    double plane_d[3];
    double cross[3][3], reci[3][3], inv[3][3];

    total_bins = 1;

    // calculate the inverse matrix of cell and the distance between cell plane
    cross[0][0] = cell[1][1]*cell[2][2] - cell[1][2]*cell[2][1];
    cross[0][1] = cell[1][2]*cell[2][0] - cell[1][0]*cell[2][2];
    cross[0][2] = cell[1][0]*cell[2][1] - cell[1][1]*cell[2][0];
    cross[1][0] = cell[2][1]*cell[0][2] - cell[2][2]*cell[0][1];
    cross[1][1] = cell[2][2]*cell[0][0] - cell[2][0]*cell[0][2];
    cross[1][2] = cell[2][0]*cell[0][1] - cell[2][1]*cell[0][0];
    cross[2][0] = cell[0][1]*cell[1][2] - cell[0][2]*cell[1][1];
    cross[2][1] = cell[0][2]*cell[1][0] - cell[0][0]*cell[1][2];
    cross[2][2] = cell[0][0]*cell[1][1] - cell[0][1]*cell[1][0];

    vol = cross[0][0]*cell[0][0] + cross[0][1]*cell[0][1] + cross[0][2]*cell[0][2];
    inv[0][0] = cross[0][0]/vol;
    inv[0][1] = cross[1][0]/vol;
    inv[0][2] = cross[2][0]/vol;
    inv[1][0] = cross[0][1]/vol;
    inv[1][1] = cross[1][1]/vol;
    inv[1][2] = cross[2][1]/vol;
    inv[2][0] = cross[0][2]/vol;
    inv[2][1] = cross[1][2]/vol;
    inv[2][2] = cross[2][2]/vol;

    int min_bins[3];
    for (int i=0; i<3; ++i) {
        double cell_dim = sqrt(cell[i][0]*cell[i][0] + cell[i][1]*cell[i][1] + cell[i][2]*cell[i][2]);
        min_bins[i] = ceil(cell_dim / cutoff);
    }

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        
        if (nbins[i] < min_bins[i]) {
            nbins[i] = min_bins[i];
        }
        
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = floor(scale[i][j] * (double) nbins[j]);
            if (bin_i[i][j] == nbins[j]) {
                bin_i[i][j]--;
            }
        }

        bin_i[i][3] = bin_i[i][0] + nbins[0]*bin_i[i][1] + nbins[0]*nbins[1]*bin_i[i][2];
        atoms_bin[bin_i[i][3]]++;
    }

    max_atoms_bin = 0;
    for (int i=0; i < total_bins; ++i) {
        if (atoms_bin[i] > max_atoms_bin)
            max_atoms_bin = atoms_bin[i];
    }

    delete[] atoms_bin;

    // # of bins in each direction
    neigh_check_bins = 1;
    for (int i=0; i < 3; ++i) {
        bin_range[i] = ceil(cutoff * nbins[i] / plane_d[i]);
        neigh_check_bins *= 2*bin_range[i];
    }
}

void find_neighbors(double** cell, double** cart, double* occupancies, double** ref_cart, 
                    double** scale, double** ref_scale, int* pbc_bools,int* atom_i, int natoms,
                    int calc_ii, double cutoff_sqr, int* bin_range, int* nbins, int** bin_i,
                    int& nneigh, double* nei_list_d, int* nei_list_i, double* nei_list_occupancy){
    
    int bin_num;
    int cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    double tmp_r2;
    double total_shift[3];
    
    for (int j=0; j < 3; ++j) {
        int temp_bin_i = ref_scale[calc_ii][j] * (double) nbins[j];
        max_bin[j] = temp_bin_i + bin_range[j];
        min_bin[j] = temp_bin_i - bin_range[j];
    }

    for (int dx=min_bin[0]; dx < max_bin[0]+1; ++dx) {
        for (int dy=min_bin[1]; dy < max_bin[1]+1; ++dy) {
            for (int dz=min_bin[2]; dz < max_bin[2]+1; ++dz) {
                pbc_bin[0] = (dx%nbins[0] + nbins[0]) % nbins[0];
                pbc_bin[1] = (dy%nbins[1] + nbins[1]) % nbins[1];
                pbc_bin[2] = (dz%nbins[2] + nbins[2]) % nbins[2];
                cell_shift[0] = (dx-pbc_bin[0]) / nbins[0];
                cell_shift[1] = (dy-pbc_bin[1]) / nbins[1];
                cell_shift[2] = (dz-pbc_bin[2]) / nbins[2];

                bin_num = pbc_bin[0] + nbins[0]*pbc_bin[1] + nbins[0]*nbins[1]*pbc_bin[2];
                for (int j=0; j < natoms; ++j) {
                    if (bin_i[j][3] != bin_num)
                        continue;

                    // same atom
                    // if (!(cell_shift[0] || cell_shift[1] || cell_shift[2]) && (i == j))
                    //     continue;

                    // take care of pbc
                    if (!pbc_bools[0] && cell_shift[0] != 0)
                        continue;

                    if (!pbc_bools[1] && cell_shift[1] != 0)
                        continue;

                    if (!pbc_bools[2] && cell_shift[2] != 0)
                        continue;

                    for (int a=0; a < 3; ++a) {
                        total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                            + cart[j][a] - ref_cart[calc_ii][a];
                    }

                    // tmp = sqrt(total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2]);
                    tmp_r2 = total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2];

                    // if (tmp < cutoff) {
                    if (tmp_r2 < cutoff_sqr) {
                        for (int a=0; a < 3; ++a)
                            nei_list_d[nneigh*4 + a] = total_shift[a];
                        nei_list_d[nneigh*4 + 3] = tmp_r2;
                        nei_list_i[nneigh*2]    = atom_i[j];
                        nei_list_i[nneigh*2 + 1] = j;
                        nei_list_occupancy[nneigh] = occupancies[j];
                        nneigh++;
                    }
                }
            }
        }
    }
}


// ##################################################
// ####             Surface Harmonics            ####
// ##################################################




extern "C" int calculate_solid_gmp_elemental_sigma_gaussian_cutoff_noderiv_opt2(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double* C1_precompute_array, double* C2_precompute_array, double* lambda_precompute_array, double* gamma_precompute_array,
                                        double** mcsh) {
    // std::cout << "here1" << std::endl;
    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;

    cutoff = 0.0;

    for (int i = 0; i < nsigmas; ++i) {
        for (int j = 0; j < natoms; ++j){
            int atom_index = atom_i[j];
            int atom_order = element_index_to_order[atom_index];
            if (cutoff < elemental_sigma_cutoffs[i][atom_order] )
                cutoff = elemental_sigma_cutoffs[i][atom_order];
        }
    }
    // cutoff = 30.0;

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }
    // std::cout << "here2" << std::endl;
    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
    int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
    double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
    int start;
    // std::cout << "here3" << std::endl;
    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        // int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        // double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);
        
        start = 0;
        // std::cout << "here4" << std::endl;
        // std::cout << "loop" << std::endl;
        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            // int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            
            SolidGMPFunctionNoderivOpt2 mcsh_function = get_solid_mcsh_function_noderiv_opt2(mcsh_order);
            

            int num_order_values = get_num_order_values(mcsh_order);
            double* desc_values = new double[num_order_values]();
            // std::cout << "here5" << std::endl;

            for (int j = 0; j < nneigh; ++j) {
                double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                
                int neigh_atom_element_index = nei_list_i[j*2];
                int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                    continue;
                double occ = nei_list_occupancy[j];
                int precompute_access_index_const = m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian;
                for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                    double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                    if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                        continue;
                    double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                    if (B == 0.0)
                        continue;
                    double C1 = C1_precompute_array[precompute_access_index_const + g];
                    double C2 = C2_precompute_array[precompute_access_index_const + g];
                    double temp = C1 * exp(C2 * r0_sqr) * occ;
                    double lambda = lambda_precompute_array[precompute_access_index_const + g]; 
                    double gamma = gamma_precompute_array[precompute_access_index_const + g];
                    mcsh_function(x0, y0, z0, r0_sqr, temp, lambda, gamma, desc_values);
                    // sum_miu += m_desc[0]*occ;
                }
            }

            for (int jj = 0; jj < num_order_values; ++jj){
                // std::cout << start+jj << std::endl;
                mcsh[ii][start+jj] = desc_values[jj];

            }
            start += num_order_values;

            // double sum_square = get_desc_value_opt2(mcsh_order, desc_values);

            // sum_square = sum_square * weight;
            // if (square != 0){
            //     mcsh[ii][m] = sum_square;
            // }
            // else {
            //     mcsh[ii][m] = sqrt(sum_square);
            // }

            delete[] desc_values;
        }
        
    }
    delete[] nei_list_d;
    delete[] nei_list_i;
    delete[] nei_list_occupancy;


    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}






void PyInit_libmcsh(void) { } // for windows
