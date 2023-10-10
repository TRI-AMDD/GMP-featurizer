#include <math.h>
#include <stdio.h>
#include "calculate_gmpordernorm.h"
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


void find_neighbors2(double** cell, double** cart, double* occupancies, double** ref_cart, 
                    double** scale, double** ref_scale, int* pbc_bools,int* atom_i, int natoms,
                    int calc_ii, double cutoff_sqr, int* bin_range, int* nbins, int** bin_i,
                    int& nneigh, double* nei_list_d, int* nei_list_i, double* nei_list_occupancy, int nneigh_counter){
    
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
                            nei_list_d[(nneigh_counter + nneigh)*4 + a] = total_shift[a];
                        nei_list_d[(nneigh_counter + nneigh)*4 + 3] = tmp_r2;
                        nei_list_i[(nneigh_counter + nneigh)*2]    = atom_i[j];
                        nei_list_i[(nneigh_counter + nneigh)*2 + 1] = j;
                        nei_list_occupancy[(nneigh_counter + nneigh)] = occupancies[j];
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


extern "C" int calculate_surface_gmpordernorm_fp_deriv_ref(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

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

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
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

    for (int ii=0; ii < cal_num; ++ii) {
        // int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;

        for (int j=0; j < 3; ++j) {
            int temp_bin_i = ref_scale[ii][j] * (double) nbins[j];
            max_bin[j] = temp_bin_i + bin_range[j];
            min_bin[j] = temp_bin_i - bin_range[j];
            // max_bin[j] = bin_i[i][j] + bin_range[j];
            // min_bin[j] = bin_i[i][j] - bin_range[j];
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


                        // take care of pbc
                        if (!pbc_bools[0] && cell_shift[0] != 0)
                            continue;

                        if (!pbc_bools[1] && cell_shift[1] != 0)
                            continue;

                        if (!pbc_bools[2] && cell_shift[2] != 0)
                            continue;

                        for (int a=0; a < 3; ++a) {
                            total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                             + cart[j][a] - ref_cart[ii][a];
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

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], inv_rs = params_d[m][5];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double sum_square = 0.0;
            //double sum_square_derivative_x = 0.0, sum_square_derivative_y = 0.0, sum_square_derivative_z = 0.0;
            // std::cout << "------------" << std::endl;
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                GMPFunction mcsh_function = get_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                // std::cout << "\t" << group_index  << "\t"<< mcsh_type << "\t" << group_coefficient<< std::endl;

                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, m_desc, deriv);
                            sum_miu += m_desc[0];
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-6){
                    mcsh[ii][m] = 0.0;
                    for (int j = 0; j < nneigh; ++j) {

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] = 0.0;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int j = 0; j < nneigh; ++j) {

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] *= (0.5 / temp);
                    }

                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    return 0;
}

extern "C" int calculate_surface_gmpordernorm_noderiv_ref(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

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

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
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

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;

        for (int j=0; j < 3; ++j) {
            int temp_bin_i = ref_scale[ii][j] * (double) nbins[j];
            max_bin[j] = temp_bin_i + bin_range[j];
            min_bin[j] = temp_bin_i - bin_range[j];
            // max_bin[j] = bin_i[i][j] + bin_range[j];
            // min_bin[j] = bin_i[i][j] - bin_range[j];
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


                        // take care of pbc
                        if (!pbc_bools[0] && cell_shift[0] != 0)
                            continue;

                        if (!pbc_bools[1] && cell_shift[1] != 0)
                            continue;

                        if (!pbc_bools[2] && cell_shift[2] != 0)
                            continue;

                        for (int a=0; a < 3; ++a) {
                            total_shift[a] = cell_shift[0]*cell[0][a] + cell_shift[1]*cell[1][a] + cell_shift[2]*cell[2][a]
                                             + cart[j][a] - ref_cart[ii][a];
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

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], inv_rs = params_d[m][5];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double sum_square = 0.0;
            // std::cout << "------------" << std::endl;
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                GMPFunctionNoderiv mcsh_function = get_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                // std::cout << "\t" << group_index  << "\t"<< mcsh_type << "\t" << group_coefficient<< std::endl;

                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;
                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, m_desc);
                            sum_miu += m_desc[0]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, inv_rs, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    return 0;
}



// ##################################################
// ####             Solid Harmonics              ####
// ##################################################


// *********   references   *********
extern "C" int calculate_solid_gmpordernorm_noderiv_ref(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];
    // std::cout <<  "test"  << std::endl;
    // std::cout << cell[0][0] << "\t" << cell[0][1]<< "\t" << cell[0][2] <<std::endl;
    // std::cout << cell[1][0] << "\t" << cell[1][1]<< "\t" << cell[1][2] <<std::endl;
    // std::cout << cell[2][0] << "\t" << cell[2][1]<< "\t" << cell[2][2] <<std::endl;
    // std::cout << ref_cart[0][0] << "\t" << ref_cart[0][1]<< "\t" << ref_cart[0][2] <<std::endl;
    // std::cout << ref_cart[1][0] << "\t" << ref_cart[1][1]<< "\t" << ref_cart[1][2] <<std::endl;

    // std::cout << ref_scale[0][0] << "\t" << ref_scale[0][1]<< "\t" << ref_scale[0][2] <<std::endl;
    // std::cout << ref_scale[1][0] << "\t" << ref_scale[1][1]<< "\t" << ref_scale[1][2] <<std::endl;
    // for (int ii=0; ii < cal_num; ++ii) {
    //     for (int j=0; j < 3; ++j) {
    //         std::cout << ref_cart[ii][j] << "\t" << ref_scale[ii][j] <<std::endl;
    //     }
    // }

    // for (int ii=0; ii < natoms; ++ii) {
    //     for (int j=0; j < 3; ++j) {
    //         std::cout << cart[ii][j] << "\t" << scale[ii][j] <<std::endl;
    //     }
    // }

    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

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

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
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

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        for (int j=0; j < 3; ++j) {
            
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            int temp_bin_i = ref_scale[ii][j] * (double) nbins[j];
            // std::cout << ii << "\t"<< j << "\t"<< ref_scale[ii][j] <<"\t"<< temp_bin_i  << std::endl;
            // std::cout << ref_scale[ii][j] << "\t"<<temp_bin_i << std::endl;
            max_bin[j] = temp_bin_i + bin_range[j];
            min_bin[j] = temp_bin_i - bin_range[j];
            // max_bin[j] = bin_i[i][j] + bin_range[j];
            // min_bin[j] = bin_i[i][j] - bin_range[j];
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
                                             + cart[j][a] - ref_cart[ii][a];
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

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double sum_square = 0.0;
            // std::cout << "------------" << std::endl;
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                // std::cout << "\t" << group_index  << "\t"<< mcsh_type << "\t" << group_coefficient<< std::endl;

                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;
                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_occ_deriv_ref(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh_docc) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];
    // std::cout <<  "test"  << std::endl;
    // std::cout << cell[0][0] << "\t" << cell[0][1]<< "\t" << cell[0][2] <<std::endl;
    // std::cout << cell[1][0] << "\t" << cell[1][1]<< "\t" << cell[1][2] <<std::endl;
    // std::cout << cell[2][0] << "\t" << cell[2][1]<< "\t" << cell[2][2] <<std::endl;
    // std::cout << ref_cart[0][0] << "\t" << ref_cart[0][1]<< "\t" << ref_cart[0][2] <<std::endl;
    // std::cout << ref_cart[1][0] << "\t" << ref_cart[1][1]<< "\t" << ref_cart[1][2] <<std::endl;

    // std::cout << ref_scale[0][0] << "\t" << ref_scale[0][1]<< "\t" << ref_scale[0][2] <<std::endl;
    // std::cout << ref_scale[1][0] << "\t" << ref_scale[1][1]<< "\t" << ref_scale[1][2] <<std::endl;
    // for (int ii=0; ii < cal_num; ++ii) {
    //     for (int j=0; j < 3; ++j) {
    //         std::cout << ref_cart[ii][j] << "\t" << ref_scale[ii][j] <<std::endl;
    //     }
    // }

    // for (int ii=0; ii < natoms; ++ii) {
    //     for (int j=0; j < 3; ++j) {
    //         std::cout << cart[ii][j] << "\t" << scale[ii][j] <<std::endl;
    //     }
    // }

    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

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

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
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

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        for (int j=0; j < 3; ++j) {
            
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            int temp_bin_i = ref_scale[ii][j] * (double) nbins[j];
            // std::cout << ii << "\t"<< j << "\t"<< ref_scale[ii][j] <<"\t"<< temp_bin_i  << std::endl;
            // std::cout << ref_scale[ii][j] << "\t"<<temp_bin_i << std::endl;
            max_bin[j] = temp_bin_i + bin_range[j];
            min_bin[j] = temp_bin_i - bin_range[j];
            // max_bin[j] = bin_i[i][j] + bin_range[j];
            // min_bin[j] = bin_i[i][j] - bin_range[j];
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
                                             + cart[j][a] - ref_cart[ii][a];
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

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double sum_square = 0.0;
            // std::cout << "------------" << std::endl;
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                // std::cout << "\t" << group_index  << "\t"<< mcsh_type << "\t" << group_coefficient<< std::endl;

                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_docc[j] = 0.0;
                    }

                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu_docc;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    
                    double dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;

            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                    }
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_fp_deriv_ref(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

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

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
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

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        for (int j=0; j < 3; ++j) {
            
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            int temp_bin_i = ref_scale[ii][j] * (double) nbins[j];
            // std::cout << ref_scale[ii][j] << "\t"<<temp_bin_i << std::endl;
            max_bin[j] = temp_bin_i + bin_range[j];
            min_bin[j] = temp_bin_i - bin_range[j];
            // max_bin[j] = bin_i[i][j] + bin_range[j];
            // min_bin[j] = bin_i[i][j] - bin_range[j];
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
                                             + cart[j][a] - ref_cart[ii][a];
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

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double sum_square = 0.0;
            //double sum_square_derivative_x = 0.0, sum_square_derivative_y = 0.0, sum_square_derivative_z = 0.0;
            // std::cout << "------------" << std::endl;
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                // std::cout << "\t" << group_index  << "\t"<< mcsh_type << "\t" << group_coefficient<< std::endl;

                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);
                        
                    }
                    
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_fp_occ_deriv_ref(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh, double** dmcsh_docc) {

    int total_bins, max_atoms_bin, bin_num, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3], cell_shift[3], max_bin[3], min_bin[3], pbc_bin[3];
    //int bin_i[natoms][4];
    double vol, tmp, tmp_r2, cutoff, cutoff_sqr;
    double plane_d[3], total_shift[3];
    double cross[3][3], reci[3][3], inv[3][3];//, powtwo[nsyms];


    // Check for not implemented mcsh type.
    for (int m=0; m < nmcsh; ++m){
        bool implemented = false;
        for (int i=0; i < NUM_IMPLEMENTED_TYPE; i++) {
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i][0] && params_i[m][1] == IMPLEMENTED_MCSH_TYPE[i][1]) {
                implemented = true;
                break;
            }
        }
        if (!implemented) return 1;
    }

    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

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

    // bin: number of repetitive cells?
    for (int i=0; i<3; ++i) {
        tmp = 0;
        for (int j=0; j<3; ++j) {
            reci[i][j] = cross[i][j]/vol;
            tmp += reci[i][j]*reci[i][j];
        }
        plane_d[i] = 1/sqrt(tmp);
        nbins[i] = ceil(plane_d[i]/cutoff);
        total_bins *= nbins[i];
    }

    int *atoms_bin = new int[total_bins];
    for (int i=0; i<total_bins; ++i)
        atoms_bin[i] = 0;

    // assign the bin index to each atom
    for (int i=0; i<natoms; ++i) {
        for (int j=0; j<3; ++j) {
            bin_i[i][j] = scale[i][j] * (double) nbins[j];
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

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // int i=cal_atoms[ii];
        // calculate neighbor atoms
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        for (int j=0; j < 3; ++j) {
            
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            int temp_bin_i = ref_scale[ii][j] * (double) nbins[j];
            // std::cout << ref_scale[ii][j] << "\t"<<temp_bin_i << std::endl;
            max_bin[j] = temp_bin_i + bin_range[j];
            min_bin[j] = temp_bin_i - bin_range[j];
            // max_bin[j] = bin_i[i][j] + bin_range[j];
            // min_bin[j] = bin_i[i][j] - bin_range[j];
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
                                             + cart[j][a] - ref_cart[ii][a];
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

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            // double weight = params_d[m][1];
            double sum_square = 0.0;
            //double sum_square_derivative_x = 0.0, sum_square_derivative_y = 0.0, sum_square_derivative_z = 0.0;
            // std::cout << "------------" << std::endl;
            // std::cout << mcsh_order << "\t" << num_groups  << std::endl;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                // std::cout << "\t" << group_index  << "\t"<< mcsh_type << "\t" << group_coefficient<< std::endl;

                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;

                        sum_dmiu_docc[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    delete [] sum_dmiu_docc;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz, dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;

                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);

                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                        
                    }
                    
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}




// *********   with preset cutoff   *********

extern "C" int calculate_solid_gmpordernorm_noderiv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh) {

    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;

    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);


    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {

        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;
                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_occ_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh_docc) {

    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;

    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_docc[j] = 0.0;
                    }

                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu_docc;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    
                    double dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;

            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                    }
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_fp_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh) {

    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;

    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);
                        
                    }
                    
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_fp_occ_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh, double** dmcsh_docc) {

    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;

    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {

        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;

                        sum_dmiu_docc[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {

                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    delete [] sum_dmiu_docc;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz, dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;

                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);

                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                        
                    }
                    
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}




// *********   with cutoff based on probe sigma   *********

extern "C" int calculate_solid_gmpordernorm_sigma_cutoff_noderiv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh) {

    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;


    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], mcsh_cutoff_sqr = params_d[m][4] * params_d[m][4];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                if (mcsh_type == 1){
                    double sum_miu = 0.0;
                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_sigma_cutoff_occ_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh_docc) {
    
    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;

    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], mcsh_cutoff_sqr = params_d[m][4] * params_d[m][4];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_docc[j] = 0.0;
                    }

                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu_docc;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    
                    double dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;

            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                    }
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_sigma_cutoff_fp_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh) {
    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;

    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], mcsh_cutoff_sqr = params_d[m][4] * params_d[m][4];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);
                        
                    }
                    
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_sigma_cutoff_fp_occ_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, int* element_index_to_order,
                                        double** mcsh, double** dmcsh, double** dmcsh_docc) {
    double cutoff, cutoff_sqr;

    // Check for not implemented mcsh type.
    if (!check_implementation(nmcsh,params_i)) return 1;

    cutoff = 0.0;
    // let cutoff equal to the maximum of Rc
    for (int m = 0; m < nmcsh; ++m) {
        if (cutoff < params_d[m][4])
            cutoff = params_d[m][4];
    }

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {

        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3], mcsh_cutoff_sqr = params_d[m][4] * params_d[m][4];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;

                        sum_dmiu_docc[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    delete [] sum_dmiu_docc;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        if (r0_sqr > mcsh_cutoff_sqr)
                            continue;
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz, dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;

                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);

                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                        
                    }
                    
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}




// *********   with cutoff based on probe and longest elemental sigma   *********

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_cutoff_noderiv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, int* element_index_to_order,
                                        double** mcsh) {

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

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {

        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                if (mcsh_type == 1){
                    double sum_miu = 0.0;
                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_cutoff_occ_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, int* element_index_to_order,
                                        double** mcsh, double** dmcsh_docc) {
    
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

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_docc[j] = 0.0;
                    }

                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu_docc;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    
                    double dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;

            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                    }
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, int* element_index_to_order,
                                        double** mcsh, double** dmcsh) {

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

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;

                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);
                        
                    }
                    
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_cutoff_fp_occ_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, int* element_index_to_order,
                                        double** mcsh, double** dmcsh, double** dmcsh_docc) {

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

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;

                        sum_dmiu_docc[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0] * occ;
                            sum_dmiu_dxj[j] += deriv[0] * occ;
                            sum_dmiu_dyj[j] += deriv[1] * occ;
                            sum_dmiu_dzj[j] += deriv[2] * occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    delete [] sum_dmiu_docc;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0] * occ;
                            sum_miu2 += miu[1] * occ;
                            sum_miu3 += miu[2] * occ;
                            sum_dmiu1_dxj[j] += deriv[0] * occ;
                            sum_dmiu1_dyj[j] += deriv[1] * occ;
                            sum_dmiu1_dzj[j] += deriv[2] * occ;
                            sum_dmiu2_dxj[j] += deriv[3] * occ;
                            sum_dmiu2_dyj[j] += deriv[4] * occ;
                            sum_dmiu2_dzj[j] += deriv[5] * occ;
                            sum_dmiu3_dxj[j] += deriv[6] * occ;
                            sum_dmiu3_dyj[j] += deriv[7] * occ;
                            sum_dmiu3_dzj[j] += deriv[8] * occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff*elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0] * occ;
                            sum_miu2 += miu[1] * occ;
                            sum_miu3 += miu[2] * occ;
                            sum_miu4 += miu[3] * occ;
                            sum_miu5 += miu[4] * occ;
                            sum_miu6 += miu[5] * occ;
                            sum_dmiu1_dxj[j] += deriv[0] * occ;
                            sum_dmiu1_dyj[j] += deriv[1] * occ;
                            sum_dmiu1_dzj[j] += deriv[2] * occ;
                            sum_dmiu2_dxj[j] += deriv[3] * occ;
                            sum_dmiu2_dyj[j] += deriv[4] * occ;
                            sum_dmiu2_dzj[j] += deriv[5] * occ;
                            sum_dmiu3_dxj[j] += deriv[6] * occ;
                            sum_dmiu3_dyj[j] += deriv[7] * occ;
                            sum_dmiu3_dzj[j] += deriv[8] * occ;
                            sum_dmiu4_dxj[j] += deriv[9] * occ;
                            sum_dmiu4_dyj[j] += deriv[10] * occ;
                            sum_dmiu4_dzj[j] += deriv[11] * occ;
                            sum_dmiu5_dxj[j] += deriv[12] * occ;
                            sum_dmiu5_dyj[j] += deriv[13] * occ;
                            sum_dmiu5_dzj[j] += deriv[14] * occ;
                            sum_dmiu6_dxj[j] += deriv[15] * occ;
                            sum_dmiu6_dyj[j] += deriv[16] * occ;
                            sum_dmiu6_dzj[j] += deriv[17] * occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz, dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;

                        // dmcsh[ii*nmcsh + m][i*3]     -= dmdx;
                        // dmcsh[ii*nmcsh + m][i*3 + 1] -= dmdy;
                        // dmcsh[ii*nmcsh + m][i*3 + 2] -= dmdz;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;

            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;

                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);

                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                        
                    }
                    
                }
            }


            // if (square != 0){
            //     mcsh[ii][m] = sum_square;
            // }
            // else {
            //     double temp = sqrt(sum_square);
            //     if (fabs(temp) < 1e-8){
            //         mcsh[ii][m] = 0.0;
            //         for (int j = 0; j < nneigh; ++j) {

            //             dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] = 0.0;
            //             dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] = 0.0;
            //             dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] = 0.0;

            //             dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] = 0.0;
            //         }
            //     }
            //     else {
            //         mcsh[ii][m] = temp;
            //         for (int j = 0; j < nneigh; ++j) {

            //             dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] *= (0.5 / temp);
            //             dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] *= (0.5 / temp);
            //             dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] *= (0.5 / temp);

            //             dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] *= (0.5 / temp);
            //         }

            //     }
            // }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}




// *********   with cutoff based on probe and each elemental sigma   *********
extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double** mcsh) {

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

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                if (mcsh_type == 1){
                    double sum_miu = 0.0;
                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        assert(neigh_atom_element_index >= 0 && neigh_atom_element_index < 120);
                        assert(neigh_atom_element_order >= 0 && neigh_atom_element_order < 120);
                        assert(j*2 < max_atoms_bin * 2 * neigh_check_bins);
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            assert(sigma_index < nmcsh && sigma_index >= 0);
                            assert(neigh_atom_element_order*max_n_gaussian+g < 119 * max_n_gaussian);
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            assert(g*2+1 < max_n_gaussian*2);
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        assert(neigh_atom_element_index >= 0 && neigh_atom_element_index < 120);
                        assert(neigh_atom_element_order >= 0 && neigh_atom_element_order < 120);
                        assert(j*2 < max_atoms_bin * 2 * neigh_check_bins);
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            assert(sigma_index < nmcsh && sigma_index >= 0);
                            assert(neigh_atom_element_order*max_n_gaussian+g < 119 * max_n_gaussian);
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        assert(neigh_atom_element_index >= 0 && neigh_atom_element_index < 120);
                        assert(neigh_atom_element_order >= 0 && neigh_atom_element_order < 120);
                        assert(j*2 < max_atoms_bin * 2 * neigh_check_bins);
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            assert(sigma_index < nmcsh && sigma_index >= 0);
                            assert(neigh_atom_element_order*max_n_gaussian+g < 119 * max_n_gaussian);
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}



// *********   with cutoff based on probe and each elemental sigma   *********
extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double* C1_precompute_array, double* C2_precompute_array, double* lambda_precompute_array, double* gamma_precompute_array,
                                        double** mcsh) {

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

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);
        // std::cout << "precompute" << std::endl;
        double* temp_precompute_array = new double[nmcsh * nneigh * max_n_gaussian];
        // for (int m = 0; m < nmcsh; ++m) {
        //     int sigma_index = params_i[m][3];
        //     for (int j = 0; j < nneigh; ++j) {
        //         int neigh_atom_element_index = nei_list_i[j*2];
        //         int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
        //         for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
        //             double r0_sqr = nei_list_d[j*4+3];
        //             double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
        //             if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
        //                 continue;
        //             double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
        //             if (B == 0.0)
        //                 continue;
        //             double C1 = C1_precompute_array[m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian + g];
        //             double C2 = C2_precompute_array[m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian + g];
                    
        //             temp_precompute_array[m * nneigh * max_n_gaussian + j * max_n_gaussian + g] = C1 * exp(C2 * r0_sqr); 
        //         }
        //     }
        // }

        for (int m = 0; m < nmcsh; ++m) {
            int sigma_index = params_i[m][3];
            int temp_precompute_access_index_m = m * nneigh * max_n_gaussian;
            int precompute_access_index_m = m * 120 * max_n_gaussian;

            for (int j = 0; j < nneigh; ++j) {
                int neigh_atom_element_index = nei_list_i[j*2];
                int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                int temp_precompute_access_index_j = temp_precompute_access_index_m + j * max_n_gaussian;
                
                for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g) {
                    double r0_sqr = nei_list_d[j*4+3];
                    double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                    double cutoff_squared = elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff;

                    if (r0_sqr > cutoff_squared)
                        continue;

                    double B = atom_gaussian[neigh_atom_element_order][g*2];
                    if (B == 0.0)
                        continue;

                    double C1 = C1_precompute_array[precompute_access_index_m + neigh_atom_element_order * max_n_gaussian + g];
                    double C2 = C2_precompute_array[precompute_access_index_m + neigh_atom_element_order * max_n_gaussian + g];
                    temp_precompute_array[temp_precompute_access_index_j + g] = C1 * exp(C2 * r0_sqr); 
                }
            }
        }

        // std::cout << "loop" << std::endl;
        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;

            // std::cout << "loop1" << std::endl;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderivOpt mcsh_function = get_solid_mcsh_function_noderiv_opt(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                if (mcsh_type == 1){
                    double sum_miu = 0.0;
                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        int temp_precompute_access_index_const = m * nneigh * max_n_gaussian + j * max_n_gaussian;
                        int precompute_access_index_const = m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian;
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            double temp = temp_precompute_array[temp_precompute_access_index_const + g];
                            double lambda = lambda_precompute_array[precompute_access_index_const + g]; 
                            double gamma = gamma_precompute_array[precompute_access_index_const + g];
                            mcsh_function(x0, y0, z0, r0_sqr, temp, lambda, gamma, m_desc);
                            sum_miu += m_desc[0]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        int temp_precompute_access_index_const = m * nneigh * max_n_gaussian + j * max_n_gaussian;
                        int precompute_access_index_const = m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian;
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            double temp = temp_precompute_array[temp_precompute_access_index_const + g];
                            double lambda = lambda_precompute_array[precompute_access_index_const + g]; 
                            double gamma = gamma_precompute_array[precompute_access_index_const + g];
                            mcsh_function(x0, y0, z0, r0_sqr, temp, lambda, gamma, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        int temp_precompute_access_index_const = m * nneigh * max_n_gaussian + j * max_n_gaussian;
                        int precompute_access_index_const = m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian;
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            double temp = temp_precompute_array[temp_precompute_access_index_const + g];
                            double lambda = lambda_precompute_array[precompute_access_index_const + g]; 
                            double gamma = gamma_precompute_array[precompute_access_index_const + g];
                            mcsh_function(x0, y0, z0, r0_sqr, temp, lambda, gamma, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                }
            }
            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
        delete[] temp_precompute_array;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}





extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt2(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double* C1_precompute_array, double* C2_precompute_array, double* lambda_precompute_array, double* gamma_precompute_array,
                                        double** mcsh) {

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

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
    int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
    double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        // int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        // double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);
        

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

            double sum_square = get_desc_value_opt2(mcsh_order, desc_values);

            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

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


extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt2_2(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double* C1_precompute_array, double* C2_precompute_array, double* lambda_precompute_array, double* gamma_precompute_array,
                                        double** mcsh) {

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

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
    int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
    double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        // int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        // double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);
        

        // std::cout << "loop" << std::endl;
        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            // int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            
            SolidGMPFunctionNoderivOpt2_2 mcsh_function = get_solid_mcsh_function_noderiv_opt2_2(mcsh_order);
            

            int num_order_values = get_num_order_values(mcsh_order);
            double* desc_values = new double[num_order_values]();
            double* temp_desc_values = new double[num_order_values]();

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
                    mcsh_function(x0, y0, z0, r0_sqr, temp, lambda, gamma, temp_desc_values);
                    for (int iii = 0; iii < num_order_values; ++iii) {
                        desc_values[iii] += temp_desc_values[iii];
                    }
                    // sum_miu += m_desc[0]*occ;
                }
            }

            double sum_square = get_desc_value_opt2(mcsh_order, desc_values);

            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                mcsh[ii][m] = sqrt(sum_square);
            }

            delete[] desc_values;
            delete[] temp_desc_values;
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


// extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt3(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
//                                         int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
//                                         int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
//                                         double* C1_precompute_array, double* C2_precompute_array, double* lambda_precompute_array, double* gamma_precompute_array,
//                                         double** mcsh) {

//     double cutoff, cutoff_sqr;

//     // Check for not implemented mcsh type.
//     if (!check_implementation(nmcsh,params_i)) return 1;

//     cutoff = 0.0;

//     for (int i = 0; i < nsigmas; ++i) {
//         for (int j = 0; j < natoms; ++j){
//             int atom_index = atom_i[j];
//             int atom_order = element_index_to_order[atom_index];
//             if (cutoff < elemental_sigma_cutoffs[i][atom_order] )
//                 cutoff = elemental_sigma_cutoffs[i][atom_order];
//         }
//     }
//     // cutoff = 30.0;

//     cutoff_sqr = cutoff * cutoff;

//     int max_atoms_bin, neigh_check_bins, nneigh;
//     int bin_range[3], nbins[3]; 
//     int **bin_i = new int*[natoms];
//     for (int i=0; i<natoms; i++) {
//         bin_i[i] = new int[4];
//     }

//     calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

//     double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
//     int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
//     double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
//     int nneigh_gaussian;
//     //for (int i=0; i < natoms; ++i) {
//     for (int ii=0; ii < cal_num; ++ii) {
//         // double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
//         // int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
//         // double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
//         nneigh = 0;
        
//         find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
//                        ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);
        
//         double* C1_array = new double[nneigh * max_n_gaussian];
//         double* C2_array = new double[nneigh * max_n_gaussian];
//         double* lambda_array = new double[nneigh * max_n_gaussian];
//         double* gamma_array = new double[nneigh * max_n_gaussian];
//         double* x0_array = new double[nneigh * max_n_gaussian];
//         double* y0_array = new double[nneigh * max_n_gaussian];
//         double* z0_array = new double[nneigh * max_n_gaussian];
//         double* r0sqr_array = new double[nneigh * max_n_gaussian];
//         double* occ_array = new double[nneigh * max_n_gaussian];

        
//         // std::cout << "loop" << std::endl;
//         for (int m = 0; m < nmcsh; ++m) {
//             int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
//             // int num_groups = get_num_groups(mcsh_order);
//             // params_d: sigma, weight, A, alpha, inv_rs
//             double A = params_d[m][2], alpha = params_d[m][3];
//             double weight = 1.0;
            
//             SolidGMPFunctionNoderivOpt3 mcsh_function = get_solid_mcsh_function_noderiv_opt3(mcsh_order);
            

//             // int num_order_values = get_num_order_values(mcsh_order);
//             // double* desc_values = new double[num_order_values]();

//             nneigh_gaussian = 0;
//             for (int j = 0; j < nneigh; ++j) {
//                 double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                
//                 int neigh_atom_element_index = nei_list_i[j*2];
//                 int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
//                 double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
//                 if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
//                     continue;
//                 double occ = nei_list_occupancy[j];
//                 int precompute_access_index_const = m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian;
//                 for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
//                     double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
//                     if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
//                         continue;
//                     double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
//                     if (B == 0.0)
//                         continue;

//                     C1_array[nneigh_gaussian] = C1_precompute_array[precompute_access_index_const + g];
//                     C2_array[nneigh_gaussian] = C2_precompute_array[precompute_access_index_const + g];
//                     lambda_array[nneigh_gaussian] = lambda_precompute_array[precompute_access_index_const + g];
//                     gamma_array[nneigh_gaussian] = gamma_precompute_array[precompute_access_index_const + g];
//                     x0_array[nneigh_gaussian] = x0;
//                     y0_array[nneigh_gaussian] = y0;
//                     z0_array[nneigh_gaussian] = z0;
//                     r0sqr_array[nneigh_gaussian] = r0_sqr;
//                     occ_array[nneigh_gaussian] = occ;
//                     nneigh_gaussian++;
//                     // double C1 = C1_precompute_array[precompute_access_index_const + g];
//                     // double C2 = C2_precompute_array[precompute_access_index_const + g];
//                     // double temp = C1 * exp(C2 * r0_sqr) * occ;
//                     // double lambda = lambda_precompute_array[precompute_access_index_const + g]; 
//                     // double gamma = gamma_precompute_array[precompute_access_index_const + g];
//                     // mcsh_function(x0, y0, z0, r0_sqr, temp, lambda, gamma, desc_values);
//                     // sum_miu += m_desc[0]*occ;
//                 }
//             }

//             // double sum_square = get_desc_value_opt2(mcsh_order, desc_values);
//             double sum_square = mcsh_function(x0_array, y0_array, z0_array, r0sqr_array, occ_array, C1_array, C2_array, lambda_array, gamma_array, nneigh_gaussian);

//             // sum_square = sum_square * weight;
//             if (square != 0){
//                 mcsh[ii][m] = sum_square;
//             }
//             else {
//                 mcsh[ii][m] = sqrt(sum_square);
//             }

//             // delete[] desc_values;
//         }

//         delete[] C1_array;
//         delete[] C2_array;
//         delete[] lambda_array;
//         delete[] gamma_array;
//         delete[] x0_array;
//         delete[] y0_array;
//         delete[] z0_array;
//         delete[] r0sqr_array;
//         delete[] occ_array;
        
//     }
//     delete[] nei_list_d;
//     delete[] nei_list_i;
//     delete[] nei_list_occupancy;


//     for (int i=0; i<natoms; i++) {
//         delete[] bin_i[i];
//     }
//     delete[] bin_i;
//     return 0;
// }



// extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt3_2_back(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
//                                         int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
//                                         int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
//                                         double* C1_precompute_array, double* C2_precompute_array, double* lambda_precompute_array, double* gamma_precompute_array,
//                                         double** mcsh) {

//     double cutoff, cutoff_sqr;

//     // Check for not implemented mcsh type.
//     if (!check_implementation(nmcsh,params_i)) return 1;

//     cutoff = 0.0;

//     for (int i = 0; i < nsigmas; ++i) {
//         for (int j = 0; j < natoms; ++j){
//             int atom_index = atom_i[j];
//             int atom_order = element_index_to_order[atom_index];
//             if (cutoff < elemental_sigma_cutoffs[i][atom_order] )
//                 cutoff = elemental_sigma_cutoffs[i][atom_order];
//         }
//     }
//     // cutoff = 30.0;

//     cutoff_sqr = cutoff * cutoff;

//     int max_atoms_bin, neigh_check_bins, nneigh;
//     int bin_range[3], nbins[3]; 
//     int **bin_i = new int*[natoms];
//     for (int i=0; i<natoms; i++) {
//         bin_i[i] = new int[4];
//     }

//     calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

//     double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
//     int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
//     double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];

//     // double* C1_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];
//     // double* C2_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];
//     // double* lambda_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];
//     // double* gamma_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];
//     // double* x0_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];
//     // double* y0_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];
//     // double* z0_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];
//     // double* r0sqr_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];
//     // double* occ_array = new double[max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num];

//     // double* C1_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     // double* C2_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     // double* lambda_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     // double* gamma_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     // double* x0_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     // double* y0_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     // double* z0_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     // double* r0sqr_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     // double* occ_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);

//     int nneigh_gaussian;

//     // std::cout << "loop" << std::endl;
//     for (int m = 0; m < nmcsh; ++m) {
//         int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
//         // int num_groups = get_num_groups(mcsh_order);
//         // params_d: sigma, weight, A, alpha, inv_rs
//         double A = params_d[m][2], alpha = params_d[m][3];
//         double weight = 1.0;
        
//         SolidGMPFunctionNoderivOpt3 mcsh_function = get_solid_mcsh_function_noderiv_opt3(mcsh_order);
//         nneigh_gaussian = 0;
//         //for (int i=0; i < natoms; ++i) {
//         for (int ii=0; ii < cal_num; ++ii) {
//             // double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
//             // int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
//             // double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
//             nneigh = 0;
            
//             find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
//                         ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

//             // int num_order_values = get_num_order_values(mcsh_order);
//             // double* desc_values = new double[num_order_values]();

            
//             for (int j = 0; j < nneigh; ++j) {
//                 double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                
//                 int neigh_atom_element_index = nei_list_i[j*2];
//                 int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
//                 double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
//                 if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
//                     continue;
//                 double occ = nei_list_occupancy[j];
//                 int precompute_access_index_const = m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian;
//                 for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
//                     double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
//                     if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
//                         continue;
//                     double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
//                     if (B == 0.0)
//                         continue;

//                     // C1_array[nneigh_gaussian] = C1_precompute_array[precompute_access_index_const + g];
//                     // C2_array[nneigh_gaussian] = C2_precompute_array[precompute_access_index_const + g];
//                     // lambda_array[nneigh_gaussian] = lambda_precompute_array[precompute_access_index_const + g];
//                     // gamma_array[nneigh_gaussian] = gamma_precompute_array[precompute_access_index_const + g];
//                     // x0_array[nneigh_gaussian] = x0;
//                     // y0_array[nneigh_gaussian] = y0;
//                     // z0_array[nneigh_gaussian] = z0;
//                     // r0sqr_array[nneigh_gaussian] = r0_sqr;
//                     // occ_array[nneigh_gaussian] = occ;
//                     // nneigh_gaussian++;
//                     // double C1 = C1_precompute_array[precompute_access_index_const + g];
//                     // double C2 = C2_precompute_array[precompute_access_index_const + g];
//                     // double temp = C1 * exp(C2 * r0_sqr) * occ;
//                     // double lambda = lambda_precompute_array[precompute_access_index_const + g]; 
//                     // double gamma = gamma_precompute_array[precompute_access_index_const + g];
//                     // mcsh_function(x0, y0, z0, r0_sqr, temp, lambda, gamma, desc_values);
//                     // sum_miu += m_desc[0]*occ;
//                 }
//             }

//             // // double sum_square = get_desc_value_opt2(mcsh_order, desc_values);
//             // double sum_square = mcsh_function(x0_array, y0_array, z0_array, r0sqr_array, occ_array, C1_array, C2_array, lambda_array, gamma_array, nneigh_gaussian);

//             // // sum_square = sum_square * weight;
//             // if (square != 0){
//             //     mcsh[ii][m] = sum_square;
//             // }
//             // else {
//             //     mcsh[ii][m] = sqrt(sum_square);
//             // }

//             // delete[] desc_values;
//         }

//         // mcsh_function(x0_array, y0_array, z0_array, r0sqr_array, occ_array, C1_array, C2_array, lambda_array, gamma_array, nneigh_gaussian);

        
        
//     }

//     // delete[] C1_array;
//     // delete[] C2_array;
//     // delete[] lambda_array;
//     // delete[] gamma_array;
//     // delete[] x0_array;
//     // delete[] y0_array;
//     // delete[] z0_array;
//     // delete[] r0sqr_array;
//     // delete[] occ_array;

//     // mkl_free(C1_array);
//     // mkl_free(C2_array);
//     // mkl_free(lambda_array);
//     // mkl_free(gamma_array);
//     // mkl_free(x0_array);
//     // mkl_free(y0_array);
//     // mkl_free(z0_array);
//     // mkl_free(r0sqr_array);
//     // mkl_free(occ_array);


//     delete[] nei_list_d;
//     delete[] nei_list_i;
//     delete[] nei_list_occupancy;


//     for (int i=0; i<natoms; i++) {
//         delete[] bin_i[i];
//     }
//     delete[] bin_i;
//     return 0;
// }


// extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_noderiv_opt3_2(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
//                                         int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
//                                         int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
//                                         double* C1_precompute_array, double* C2_precompute_array, double* lambda_precompute_array, double* gamma_precompute_array,
//                                         double** mcsh) {

//     double cutoff, cutoff_sqr;

//     // Check for not implemented mcsh type.
//     if (!check_implementation(nmcsh,params_i)) return 1;

//     cutoff = 0.0;

//     for (int i = 0; i < nsigmas; ++i) {
//         for (int j = 0; j < natoms; ++j){
//             int atom_index = atom_i[j];
//             int atom_order = element_index_to_order[atom_index];
//             if (cutoff < elemental_sigma_cutoffs[i][atom_order] )
//                 cutoff = elemental_sigma_cutoffs[i][atom_order];
//         }
//     }
//     // cutoff = 30.0;

//     cutoff_sqr = cutoff * cutoff;

//     int max_atoms_bin, neigh_check_bins, nneigh;
//     int bin_range[3], nbins[3]; 
//     int **bin_i = new int*[natoms];
//     for (int i=0; i<natoms; i++) {
//         bin_i[i] = new int[4];
//     }

//     calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

//     double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins * cal_num];
//     int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins * cal_num];
//     double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins * cal_num];

//     int*    nneigh_list = new int[cal_num];
//     int*    nneigh_counter_list = new int[cal_num];
//     int nneigh_counter = 0;

//     for (int ii=0; ii < cal_num; ++ii) {
//         // double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
//         // int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
//         // double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
//         nneigh = 0;
//         nneigh_counter_list[ii] = nneigh_counter;
//         find_neighbors2(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
//                        ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy, nneigh_counter);

//         nneigh_list[ii] = nneigh;
//         nneigh_counter += nneigh;

//     }

//     double* C1_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     double* C2_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     double* lambda_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     double* gamma_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     double* x0_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     double* y0_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     double* z0_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     double* r0sqr_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     double* occ_array = (double*)mkl_malloc(max_atoms_bin * neigh_check_bins * max_n_gaussian * cal_num * sizeof(double), 64);
//     int nneigh_gaussian;

//     // std::cout << "loop" << std::endl;
//     for (int m = 0; m < nmcsh; ++m) {
//         int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
//         // int num_groups = get_num_groups(mcsh_order);
//         // params_d: sigma, weight, A, alpha, inv_rs
//         double A = params_d[m][2], alpha = params_d[m][3];
//         double weight = 1.0;
        
//         SolidGMPFunctionNoderivOpt3 mcsh_function = get_solid_mcsh_function_noderiv_opt3(mcsh_order);
//         //for (int i=0; i < natoms; ++i) {
//         for (int ii=0; ii < cal_num; ++ii) {
//             // double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
//             // int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
//             // double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
//             nneigh = nneigh_list[ii];
//             nneigh_counter = nneigh_counter_list[ii];
        
//             // find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
//             //                ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);


        
        
            

//             // int num_order_values = get_num_order_values(mcsh_order);
//             // double* desc_values = new double[num_order_values]();

//             nneigh_gaussian = 0;
//             for (int j = 0; j < nneigh; ++j) {
//                 double x0 = nei_list_d[(nneigh_counter+j)*4], y0 = nei_list_d[(nneigh_counter+j)*4+1], z0 = nei_list_d[(nneigh_counter+j)*4+2], r0_sqr = nei_list_d[(nneigh_counter+j)*4+3];
                
//                 int neigh_atom_element_index = nei_list_i[(nneigh_counter+j)*2];
//                 int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
//                 double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
//                 if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
//                     continue;
//                 double occ = nei_list_occupancy[j];
//                 int precompute_access_index_const = m * 120 * max_n_gaussian + neigh_atom_element_order * max_n_gaussian;
//                 for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
//                     double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
//                     if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
//                         continue;
//                     double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
//                     if (B == 0.0)
//                         continue;

//                     C1_array[nneigh_gaussian] = C1_precompute_array[precompute_access_index_const + g];
//                     C2_array[nneigh_gaussian] = C2_precompute_array[precompute_access_index_const + g];
//                     lambda_array[nneigh_gaussian] = lambda_precompute_array[precompute_access_index_const + g];
//                     gamma_array[nneigh_gaussian] = gamma_precompute_array[precompute_access_index_const + g];
//                     x0_array[nneigh_gaussian] = x0;
//                     y0_array[nneigh_gaussian] = y0;
//                     z0_array[nneigh_gaussian] = z0;
//                     r0sqr_array[nneigh_gaussian] = r0_sqr;
//                     occ_array[nneigh_gaussian] = occ;
//                     nneigh_gaussian++;

//                 }
//             }

//             // double sum_square = get_desc_value_opt2(mcsh_order, desc_values);
//             // double sum_square = mcsh_function(x0_array, y0_array, z0_array, r0sqr_array, occ_array, C1_array, C2_array, lambda_array, gamma_array, nneigh_gaussian);

//             // // sum_square = sum_square * weight;
//             // if (square != 0){
//             //     mcsh[ii][m] = sum_square;
//             // }
//             // else {
//             //     mcsh[ii][m] = sqrt(sum_square);
//             // }

//             // delete[] desc_values;
//         }
//         mcsh_function(x0_array, y0_array, z0_array, r0sqr_array, occ_array, C1_array, C2_array, lambda_array, gamma_array, nneigh_gaussian);

//         // delete[] C1_array;
//         // delete[] C2_array;
//         // delete[] lambda_array;
//         // delete[] gamma_array;
//         // delete[] x0_array;
//         // delete[] y0_array;
//         // delete[] z0_array;
//         // delete[] r0sqr_array;
//         // delete[] occ_array;
        
//     }

//     mkl_free(C1_array);
//     mkl_free(C2_array);
//     mkl_free(lambda_array);
//     mkl_free(gamma_array);
//     mkl_free(x0_array);
//     mkl_free(y0_array);
//     mkl_free(z0_array);
//     mkl_free(r0sqr_array);
//     mkl_free(occ_array);

//     delete[] nei_list_d;
//     delete[] nei_list_i;
//     delete[] nei_list_occupancy;


//     for (int i=0; i<natoms; i++) {
//         delete[] bin_i[i];
//     }
//     delete[] bin_i;
//     return 0;
// }

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_occ_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double** mcsh, double** dmcsh_docc ) {
    
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

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunctionNoderiv mcsh_function = get_solid_mcsh_function_noderiv(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_docc[j] = 0.0;
                    }

                    double m_desc[1];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc);
                            sum_miu += m_desc[0]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu_docc;
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff = elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g];
                            if (r0_sqr > (elemental_sigma_gausisan_cutoff * elemental_sigma_gausisan_cutoff))
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0.0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    
                    double dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;

            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                    }
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}


extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_occ_deriv_opt2_2(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double* C1_precompute_array, double* C2_precompute_array, double* lambda_precompute_array, double* gamma_precompute_array,
                                        double** mcsh, double** dmcsh_docc) {

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

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
    int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
    double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        // double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        // int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        // double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);
        

        // std::cout << "loop" << std::endl;
        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            // int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            
            SolidGMPFunctionNoderivOpt2_2 mcsh_function = get_solid_mcsh_function_noderiv_opt2_2(mcsh_order);
            

            int num_order_values = get_num_order_values(mcsh_order);
            double* desc_values = new double[num_order_values]();
            double* temp_desc_values = new double[num_order_values]();
            double* temp_occ_sum_values = new double[nneigh * num_order_values]();

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
                    double temp = C1 * exp(C2 * r0_sqr) ;
                    double lambda = lambda_precompute_array[precompute_access_index_const + g]; 
                    double gamma = gamma_precompute_array[precompute_access_index_const + g];
                    mcsh_function(x0, y0, z0, r0_sqr, temp, lambda, gamma, temp_desc_values);
                    for (int iii = 0; iii < num_order_values; ++iii) {
                        desc_values[iii] += temp_desc_values[iii] * occ;
                        temp_occ_sum_values[j * num_order_values + iii] += temp_desc_values[iii];
                    }
                    // sum_miu += m_desc[0]*occ;
                }
            }

            int desc_idx = ii*nmcsh + m;
            double* temp_occ_sum_values2 = new double[num_order_values]();
            for (int j = 0; j < nneigh; ++j) {
                for (int iii = 0; iii < num_order_values; ++iii) {
                    temp_occ_sum_values2[iii] = temp_occ_sum_values[j * num_order_values + iii];
                }
                int atom_idx = nei_list_i[j*2 + 1];
                dmcsh_docc[desc_idx][atom_idx] = 2.0 * get_docc_value_opt2(mcsh_order, desc_values, temp_occ_sum_values2);
            }
            delete[] temp_occ_sum_values2;


            double sum_square = get_desc_value_opt2(mcsh_order, desc_values);

            

            // sum_square = sum_square * weight;
            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {
                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                    }
                }
            }


            delete[] desc_values;
            delete[] temp_desc_values;
            delete[] temp_occ_sum_values;
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


extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double** mcsh, double** dmcsh) {

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

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);
    
    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {

        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);

                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff_sqr = pow(elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g], 2);
                            if (r0_sqr > elemental_sigma_gausisan_cutoff_sqr)
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff_sqr = pow(elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g], 2);
                            if (r0_sqr > elemental_sigma_gausisan_cutoff_sqr)
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff_sqr = pow(elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g], 2);
                            if (r0_sqr > elemental_sigma_gausisan_cutoff_sqr)
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;
                }
            }
            // sum_square = sum_square * weight;

            // for (int j = 0; j < nneigh; ++j) {
            //     if (fabs(dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3]) > 100000){
            //         std::cout << dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3]  <<std::endl;
            //         std::cout << "------------" << std::endl;
            //     }               
            // }

            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {

                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3] = 0.0;
                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3 + 1] = 0.0;
                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3 + 2] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    // std::cout << ii << "\t" << m <<std::endl;
                    // std::cout << "------------" << std::endl;
                    for (int jj = 0; jj < natoms; ++jj) {
                        // if (fabs(dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3]) > 100000){
                            // std::cout << dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3] << "\t" << temp << "\t" << (0.5 / temp) * dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3] <<std::endl;
                            // std::cout << "------------" << std::endl;
                        // }
                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3] = (0.5 / temp) * dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3];
                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3 + 1] = (0.5 / temp) * dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3+1];
                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3 + 2] = (0.5 / temp) * dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3+2];
                        
                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3] = dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3];
                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3 + 1] = dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3+1];
                        // dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3 + 2] = dmcsh[ii*nmcsh + m][nei_list_i[jj*2 + 1]*3+2];


                        // std::cout << ii*nmcsh + m << "\t" << jj*3 << "\t" << jj*3+1 << "\t" << jj*3+2 <<std::endl;
                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);
                        
                    }
                    
                }
            }

            // for (int j = 0; j < nneigh; ++j) {
            //     if (fabs(dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3]) > 100000){
            //         std::cout << dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3]  <<std::endl;
            //         std::cout << "2222------------" << std::endl;
            //     }               
            // }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}

extern "C" int calculate_solid_gmpordernorm_elemental_sigma_gaussian_cutoff_fp_occ_deriv(double** cell, double** cart, double* occupancies, double** ref_cart, double** scale, double** ref_scale, int* pbc_bools,
                                        int* atom_i, int natoms, /*int* cal_atoms,*/ int cal_num, int nsigmas, int max_n_gaussian,
                                        int** params_i, double** params_d, int nmcsh, double** atom_gaussian, int* ngaussians, double** elemental_sigma_cutoffs, double** elemental_sigma_gaussian_cutoffs, int* element_index_to_order,
                                        double** mcsh, double** dmcsh, double** dmcsh_docc) {

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

    cutoff_sqr = cutoff * cutoff;

    int max_atoms_bin, neigh_check_bins, nneigh;
    int bin_range[3], nbins[3]; 
    int **bin_i = new int*[natoms];
    for (int i=0; i<natoms; i++) {
        bin_i[i] = new int[4];
    }

    calculate_bin_ranges(cell, scale, natoms, cutoff, max_atoms_bin, neigh_check_bins, bin_i, bin_range, nbins);

    //for (int i=0; i < natoms; ++i) {
    for (int ii=0; ii < cal_num; ++ii) {
        
        double* nei_list_d = new double[max_atoms_bin * 4 * neigh_check_bins];
        int*    nei_list_i = new int[max_atoms_bin * 2 * neigh_check_bins];
        double* nei_list_occupancy = new double[max_atoms_bin * neigh_check_bins];
        nneigh = 0;
        
        find_neighbors(cell, cart, occupancies, ref_cart, scale, ref_scale, pbc_bools, atom_i, natoms,
                       ii, cutoff_sqr, bin_range, nbins, bin_i, nneigh, nei_list_d, nei_list_i, nei_list_occupancy);

        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0], square = params_i[m][1], sigma_index = params_i[m][3];
            int num_groups = get_num_groups(mcsh_order);
            // params_d: sigma, weight, A, alpha, cutoff, inv_rs
            double A = params_d[m][2], alpha = params_d[m][3];
            double weight = 1.0;
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                SolidGMPFunction mcsh_function = get_solid_mcsh_function(mcsh_order, group_index);
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                
                if (mcsh_type == 1){
                    double sum_miu = 0.0;

                    double* sum_dmiu_dxj = new double[nneigh];
                    double* sum_dmiu_dyj = new double[nneigh];
                    double* sum_dmiu_dzj = new double[nneigh];
                    double* sum_dmiu_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu_dxj[j] = 0.0;
                        sum_dmiu_dyj[j] = 0.0;
                        sum_dmiu_dzj[j] = 0.0;

                        sum_dmiu_docc[j] = 0.0;
                    }
                    double m_desc[1], deriv[3];

                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff_sqr = pow(elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g], 2);
                            if (r0_sqr > elemental_sigma_gausisan_cutoff_sqr)
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, m_desc, deriv);
                            sum_miu += m_desc[0]*occ;
                            sum_dmiu_dxj[j] += deriv[0]*occ;
                            sum_dmiu_dyj[j] += deriv[1]*occ;
                            sum_dmiu_dzj[j] += deriv[2]*occ;

                            sum_dmiu_docc[j] += m_desc[0];
                        }
                    }
                    sum_square += group_coefficient * sum_miu * sum_miu;

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu * sum_dmiu_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu_dxj;
                    delete [] sum_dmiu_dyj;
                    delete [] sum_dmiu_dzj;
                    delete [] sum_dmiu_docc;
                    
                }

                if (mcsh_type == 2){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                    }

                    double miu[3], deriv[9];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff_sqr = pow(elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g], 2);
                            if (r0_sqr > elemental_sigma_gausisan_cutoff_sqr)
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);

                    double dmdx, dmdy, dmdz, dmdocc;
                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0;
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0;
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] + sum_miu3 * sum_dmiu3_docc[j]) * group_coefficient * 2.0;
                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }

                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    
                }

                if (mcsh_type == 3){
                    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;

                    double* sum_dmiu1_dxj = new double[nneigh];
                    double* sum_dmiu2_dxj = new double[nneigh];
                    double* sum_dmiu3_dxj = new double[nneigh];
                    double* sum_dmiu4_dxj = new double[nneigh];
                    double* sum_dmiu5_dxj = new double[nneigh];
                    double* sum_dmiu6_dxj = new double[nneigh];
                    double* sum_dmiu1_dyj = new double[nneigh];
                    double* sum_dmiu2_dyj = new double[nneigh];
                    double* sum_dmiu3_dyj = new double[nneigh];
                    double* sum_dmiu4_dyj = new double[nneigh];
                    double* sum_dmiu5_dyj = new double[nneigh];
                    double* sum_dmiu6_dyj = new double[nneigh];
                    double* sum_dmiu1_dzj = new double[nneigh];
                    double* sum_dmiu2_dzj = new double[nneigh];
                    double* sum_dmiu3_dzj = new double[nneigh];
                    double* sum_dmiu4_dzj = new double[nneigh];
                    double* sum_dmiu5_dzj = new double[nneigh];
                    double* sum_dmiu6_dzj = new double[nneigh];

                    double* sum_dmiu1_docc = new double[nneigh];
                    double* sum_dmiu2_docc = new double[nneigh];
                    double* sum_dmiu3_docc = new double[nneigh];
                    double* sum_dmiu4_docc = new double[nneigh];
                    double* sum_dmiu5_docc = new double[nneigh];
                    double* sum_dmiu6_docc = new double[nneigh];
                    for (int j=0; j<nneigh; j++) {
                        sum_dmiu1_dxj[j] = 0.0;
                        sum_dmiu2_dxj[j] = 0.0;
                        sum_dmiu3_dxj[j] = 0.0;
                        sum_dmiu4_dxj[j] = 0.0;
                        sum_dmiu5_dxj[j] = 0.0;
                        sum_dmiu6_dxj[j] = 0.0;
                        sum_dmiu1_dyj[j] = 0.0;
                        sum_dmiu2_dyj[j] = 0.0;
                        sum_dmiu3_dyj[j] = 0.0;
                        sum_dmiu4_dyj[j] = 0.0;
                        sum_dmiu5_dyj[j] = 0.0;
                        sum_dmiu6_dyj[j] = 0.0;
                        sum_dmiu1_dzj[j] = 0.0;
                        sum_dmiu2_dzj[j] = 0.0;
                        sum_dmiu3_dzj[j] = 0.0;
                        sum_dmiu4_dzj[j] = 0.0;
                        sum_dmiu5_dzj[j] = 0.0;
                        sum_dmiu6_dzj[j] = 0.0;

                        sum_dmiu1_docc[j] = 0.0;
                        sum_dmiu2_docc[j] = 0.0;
                        sum_dmiu3_docc[j] = 0.0;
                        sum_dmiu4_docc[j] = 0.0;
                        sum_dmiu5_docc[j] = 0.0;
                        sum_dmiu6_docc[j] = 0.0;
                    }

                    double miu[6], deriv[18];
                    for (int j = 0; j < nneigh; ++j) {
                        double x0 = nei_list_d[j*4], y0 = nei_list_d[j*4+1], z0 = nei_list_d[j*4+2], r0_sqr = nei_list_d[j*4+3];
                        int neigh_atom_element_index = nei_list_i[j*2];
                        int neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
                        double elemental_sigma_cutoff = elemental_sigma_cutoffs[sigma_index][neigh_atom_element_order];
                        if (r0_sqr > (elemental_sigma_cutoff * elemental_sigma_cutoff))
                            continue;
                        double occ = nei_list_occupancy[j];
                        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
                            double elemental_sigma_gausisan_cutoff_sqr = pow(elemental_sigma_gaussian_cutoffs[sigma_index][neigh_atom_element_order*max_n_gaussian+g], 2);
                            if (r0_sqr > elemental_sigma_gausisan_cutoff_sqr)
                                continue;
                            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
                            if (B == 0)
                                continue;
                            mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta, miu, deriv);
                            // miu: miu_1, miu_2, miu_3
                            // deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]*occ;
                            sum_miu2 += miu[1]*occ;
                            sum_miu3 += miu[2]*occ;
                            sum_miu4 += miu[3]*occ;
                            sum_miu5 += miu[4]*occ;
                            sum_miu6 += miu[5]*occ;
                            sum_dmiu1_dxj[j] += deriv[0]*occ;
                            sum_dmiu1_dyj[j] += deriv[1]*occ;
                            sum_dmiu1_dzj[j] += deriv[2]*occ;
                            sum_dmiu2_dxj[j] += deriv[3]*occ;
                            sum_dmiu2_dyj[j] += deriv[4]*occ;
                            sum_dmiu2_dzj[j] += deriv[5]*occ;
                            sum_dmiu3_dxj[j] += deriv[6]*occ;
                            sum_dmiu3_dyj[j] += deriv[7]*occ;
                            sum_dmiu3_dzj[j] += deriv[8]*occ;
                            sum_dmiu4_dxj[j] += deriv[9]*occ;
                            sum_dmiu4_dyj[j] += deriv[10]*occ;
                            sum_dmiu4_dzj[j] += deriv[11]*occ;
                            sum_dmiu5_dxj[j] += deriv[12]*occ;
                            sum_dmiu5_dyj[j] += deriv[13]*occ;
                            sum_dmiu5_dzj[j] += deriv[14]*occ;
                            sum_dmiu6_dxj[j] += deriv[15]*occ;
                            sum_dmiu6_dyj[j] += deriv[16]*occ;
                            sum_dmiu6_dzj[j] += deriv[17]*occ;

                            sum_dmiu1_docc[j] += miu[0];
                            sum_dmiu2_docc[j] += miu[1];
                            sum_dmiu3_docc[j] += miu[2];
                            sum_dmiu4_docc[j] += miu[3];
                            sum_dmiu5_docc[j] += miu[4];
                            sum_dmiu6_docc[j] += miu[5];
                        }
                    }
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
                    double dmdx, dmdy, dmdz, dmdocc;

                    for (int j = 0; j < nneigh; ++j) {
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0;

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0;

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0;

                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3] += dmdx;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 1] += dmdy;
                        dmcsh[ii*nmcsh + m][nei_list_i[j*2 + 1]*3 + 2] += dmdz;

                        dmdocc = (sum_miu1 * sum_dmiu1_docc[j] + sum_miu2 * sum_dmiu2_docc[j] +
                                  sum_miu3 * sum_dmiu3_docc[j] + sum_miu4 * sum_dmiu4_docc[j] +
                                  sum_miu5 * sum_dmiu5_docc[j] + sum_miu6 * sum_dmiu6_docc[j]) * group_coefficient * 2.0;

                        dmcsh_docc[ii*nmcsh + m][nei_list_i[j*2 + 1]] += dmdocc;
                    }
                    
                    delete [] sum_dmiu1_dxj;
                    delete [] sum_dmiu2_dxj;
                    delete [] sum_dmiu3_dxj;
                    delete [] sum_dmiu4_dxj;
                    delete [] sum_dmiu5_dxj;
                    delete [] sum_dmiu6_dxj;
                    delete [] sum_dmiu1_dyj;
                    delete [] sum_dmiu2_dyj;
                    delete [] sum_dmiu3_dyj;
                    delete [] sum_dmiu4_dyj;
                    delete [] sum_dmiu5_dyj;
                    delete [] sum_dmiu6_dyj;
                    delete [] sum_dmiu1_dzj;
                    delete [] sum_dmiu2_dzj;
                    delete [] sum_dmiu3_dzj;
                    delete [] sum_dmiu4_dzj;
                    delete [] sum_dmiu5_dzj;
                    delete [] sum_dmiu6_dzj;

                    delete [] sum_dmiu1_docc;
                    delete [] sum_dmiu2_docc;
                    delete [] sum_dmiu3_docc;
                    delete [] sum_dmiu4_docc;
                    delete [] sum_dmiu5_docc;
                    delete [] sum_dmiu6_docc;
                }
            }
            // sum_square = sum_square * weight;

            if (square != 0){
                mcsh[ii][m] = sum_square;
            }
            else {
                double temp = sqrt(sum_square);
                if (fabs(temp) < 1e-15){
                    mcsh[ii][m] = 0.0;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 1] = 0.0;
                        dmcsh[ii*nmcsh + m][jj*3 + 2] = 0.0;

                        dmcsh_docc[ii*nmcsh + m][jj] = 0.0;
                    }
                }
                else {
                    mcsh[ii][m] = temp;
                    for (int jj = 0; jj < natoms; ++jj) {

                        dmcsh[ii*nmcsh + m][jj*3] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 1] *= (0.5 / temp);
                        dmcsh[ii*nmcsh + m][jj*3 + 2] *= (0.5 / temp);

                        dmcsh_docc[ii*nmcsh + m][jj] *= (0.5 / temp);
                        
                    }
                    
                }
            }

        }
        delete[] nei_list_d;
        delete[] nei_list_i;
        delete[] nei_list_occupancy;
    }

    for (int i=0; i<natoms; i++) {
        delete[] bin_i[i];
    }
    delete[] bin_i;
    return 0;
}



void PyInit_libmcsh(void) { } // for windows
