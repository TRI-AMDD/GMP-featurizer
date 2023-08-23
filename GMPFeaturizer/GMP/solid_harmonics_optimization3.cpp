
#include "mkl.h"
#include "solid_harmonics_optimization3.h"

SolidGMPFunctionNoderivOpt3 get_solid_mcsh_function_noderiv_opt3(int mcsh_order)
{
    SolidGMPFunctionNoderivOpt3 result;
    if (mcsh_order == -1) {
        result = calc_solid_MCSH_n1_noderiv_opt3;
    } else if (mcsh_order == 0) {
        result = calc_solid_MCSH_0_noderiv_opt3;
    } else if (mcsh_order == 1) {
        result = calc_solid_MCSH_1_noderiv_opt3;
    } else if (mcsh_order == 2) {
        result = calc_solid_MCSH_2_noderiv_opt3;
    } else if (mcsh_order == 3) {
        result = calc_solid_MCSH_3_noderiv_opt3;
    // } else if (mcsh_order == 4) {
    //     result = calc_solid_MCSH_4_noderiv_opt3;
    // } else if (mcsh_order == 5) {
    //     result = calc_solid_MCSH_5_noderiv_opt3;
    // } else if (mcsh_order == 6) {
    //     result = calc_solid_MCSH_6_noderiv_opt3;
    // } else if (mcsh_order == 7) {
    //     result = calc_solid_MCSH_7_noderiv_opt3;
    // } else if (mcsh_order == 8) {
    //     result = calc_solid_MCSH_8_noderiv_opt3;
    // } else if (mcsh_order == 9) {
    //     result = calc_solid_MCSH_9_noderiv_opt3;
    }
    return result;
}







// C1 * exp(C2 * r0_sqr) * occ;

void calc_temp_opt3(double* C1, double* C2, double* r0_sqr, double* occ, double* result, int n){
    vdMul(n, C2, r0_sqr, result);
    vdExp(n, result, result);
    vdMul(n, C1, result, result);
    vdMul(n, occ, result, result);
}


// inline double P2_opt2(double lambda_x0_2, double inv_gamma){
//     return (0.5 * inv_gamma) + lambda_x0_2;
// }
void P2_vec_opt3(double* lambda_x0_2, double* inv_gamma, double* result, int n) {
    // Copy lambda_x0_2 to result
    cblas_dcopy(n, lambda_x0_2, 1, result, 1);

    // result = 0.5 * inv_gamma + result
    cblas_daxpy(n, 0.5, inv_gamma, 1, result, 1);
}


// inline double P3_opt2(double lambda_x0, double lambda_x0_3, double inv_gamma){
//     return (1.5 * lambda_x0 * inv_gamma) + lambda_x0_3;
// }

void P3_vec_opt3(double* lambda_x0, double* lambda_x0_3, double* inv_gamma, double* result, int n) {
    // Compute lambda_x0[i] * inv_gamma[i] for each i, and store directly in result
    vdMul(n, lambda_x0, inv_gamma, result);

    // result = result + 1.5 * lambda_x0_3
    cblas_daxpy(n, 1.5, lambda_x0_3, 1, result, 1);
}

// inline double P4_opt2(double lambda_x0_2, double lambda_x0_4, double inv_gamma, double inv_gamma_2){
//     return (0.75 * inv_gamma_2) + (3.0 * lambda_x0_2 * inv_gamma) + lambda_x0_4;
// }

void P4_vec_opt3(double* lambda_x0_2, double* lambda_x0_4, double* inv_gamma, double* inv_gamma_2, double* result, int n) {
    // // Compute 0.75 * inv_gamma_2[i] for each i, and store in result
    // vdMulC(n, inv_gamma_2, 0.75, result);

    // // Compute 3.0 * lambda_x0_2[i] * inv_gamma[i] for each i, and store in a temporary array
    // double* temp = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
    // vdMul(n, lambda_x0_2, inv_gamma, temp);   // temp[i] = lambda_x0_2[i] * inv_gamma[i]
    // vdMulC(n, temp, 3.0, temp);               // temp[i] = 3.0 * temp[i]

    // // Add temp[i] to result[i]
    // vdAdd(n, result, temp, result);

    // // Add lambda_x0_4[i] to each result[i]
    // vdAdd(n, lambda_x0_4, result, result);

    // // Free the temporary array
    // mkl_free(temp);

    // Compute lambda_x0_2[i] * inv_gamma[i] and store directly in result
    vdMul(n, lambda_x0_2, inv_gamma, result);
    
    // Scale result by 3.0
    cblas_dscal(n, 3.0, result, 1);
    
    // Add 0.75 * inv_gamma_2[i] to result[i]
    cblas_daxpy(n, 0.75, inv_gamma_2, 1, result, 1);

    // Add lambda_x0_4[i] to result[i]
    cblas_daxpy(n, 1.0, lambda_x0_4, 1, result, 1);
}

// // inline double P5_opt2(double lambda_x0, double lambda_x0_3, double lambda_x0_5, double inv_gamma, double inv_gamma_2){
// //     return ((15.0 * lambda_x0) * (0.25 * inv_gamma_2)) + (5.0 * lambda_x0_3 * inv_gamma) + lambda_x0_5;
// // }

// void P5_vec_opt3(double* lambda_x0, double* lambda_x0_3, double* lambda_x0_5, double* inv_gamma, double* inv_gamma_2, double* result, int n) {
//     // Compute 15.0 * lambda_x0[i] * 0.25 * inv_gamma_2[i] for each i, and store in result
//     vdMul(n, lambda_x0, inv_gamma_2, result);  // result[i] = lambda_x0[i] * inv_gamma_2[i]
//     vdMulC(n, result, 3.75, result);    // result[i] = 3.75 * result[i]

//     // Compute 5.0 * lambda_x0_3[i] * inv_gamma[i] for each i, and store in a temporary array
//     double* temp = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_3, inv_gamma, temp);    // temp[i] = lambda_x0_3[i] * inv_gamma[i]
//     vdMulC(n, temp, 5.0, temp);                // temp[i] = 5.0 * temp[i]

//     // Add temp[i] to result[i]
//     vdAdd(n, result, temp, result);

//     // Add lambda_x0_5[i] to each result[i]
//     vdAdd(n, lambda_x0_5, result, result);

//     // Free the temporary array
//     mkl_free(temp);
// }

// // inline double P6_opt2(double lambda_x0_2, double lambda_x0_4, double lambda_x0_6, double inv_gamma, double inv_gamma_2, double inv_gamma_3){
// //     return (1.875 * inv_gamma_3) + (11.25 * lambda_x0_2 * inv_gamma_2) + (7.5 * lambda_x0_4 *inv_gamma) + lambda_x0_6;
// // }

// void P6_vec_opt3(double* lambda_x0_2, double* lambda_x0_4, double* lambda_x0_6, double* inv_gamma, double* inv_gamma_2, double* inv_gamma_3, double* result, int n) {
//     // Compute 1.875 * inv_gamma_3[i] for each i, and store in result
//     vdMulC(n, inv_gamma_3, 1.875, result);

//     // Compute 11.25 * lambda_x0_2[i] * inv_gamma_2[i] for each i, and store in a temporary array
//     double* temp1 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_2, inv_gamma_2, temp1);    // temp1[i] = lambda_x0_2[i] * inv_gamma_2[i]
//     vdMulC(n, temp1, 11.25, temp1);               // temp1[i] = 11.25 * temp1[i]

//     // Compute 7.5 * lambda_x0_4[i] * inv_gamma[i] for each i, and store in another temporary array
//     double* temp2 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_4, inv_gamma, temp2);      // temp2[i] = lambda_x0_4[i] * inv_gamma[i]
//     vdMulC(n, temp2, 7.5, temp2);                 // temp2[i] = 7.5 * temp2[i]

//     // Add temp1[i] and temp2[i] to result[i]
//     vdAdd(n, result, temp1, result);
//     vdAdd(n, result, temp2, result);

//     // Add lambda_x0_6[i] to each result[i]
//     vdAdd(n, lambda_x0_6, result, result);

//     // Free the temporary arrays
//     mkl_free(temp1);
//     mkl_free(temp2);
// }

// // inline double P7_opt2(double lambda_x0, double lambda_x0_3, double lambda_x0_5, double lambda_x0_7, double inv_gamma, double inv_gamma_2, double inv_gamma_3){
// //     double term1 = 13.125 * lambda_x0 * inv_gamma_3;
// //     double term2 = 26.25 * lambda_x0_3 * inv_gamma_2;
// //     double term3 = 10.5 * lambda_x0_5 * inv_gamma;
// //     return term1 + term2 + term3 + lambda_x0_7;
// // }

// void P7_vec_opt3(double* lambda_x0, double* lambda_x0_3, double* lambda_x0_5, double* lambda_x0_7, double* inv_gamma, double* inv_gamma_2, double* inv_gamma_3, double* result, int n) {
//     // Compute 13.125 * lambda_x0[i] * inv_gamma_3[i] for each i, and store in result
//     vdMul(n, lambda_x0, inv_gamma_3, result);  // result[i] = lambda_x0[i] * inv_gamma_3[i]
//     vdMulC(n, result, 13.125, result);         // result[i] = 13.125 * result[i]

//     // Compute 26.25 * lambda_x0_3[i] * inv_gamma_2[i] for each i, and store in a temporary array
//     double* temp1 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_3, inv_gamma_2, temp1); // temp1[i] = lambda_x0_3[i] * inv_gamma_2[i]
//     vdMulC(n, temp1, 26.25, temp1);            // temp1[i] = 26.25 * temp1[i]

//     // Compute 10.5 * lambda_x0_5[i] * inv_gamma[i] for each i, and store in another temporary array
//     double* temp2 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_5, inv_gamma, temp2);   // temp2[i] = lambda_x0_5[i] * inv_gamma[i]
//     vdMulC(n, temp2, 10.5, temp2);              // temp2[i] = 10.5 * temp2[i]

//     // Add temp1[i] and temp2[i] to result[i]
//     vdAdd(n, result, temp1, result);
//     vdAdd(n, result, temp2, result);

//     // Add lambda_x0_7[i] to each result[i]
//     vdAdd(n, lambda_x0_7, result, result);

//     // Free the temporary arrays
//     mkl_free(temp1);
//     mkl_free(temp2);
// }

// // inline double P8_opt2(double lambda_x0_2, double lambda_x0_4, double lambda_x0_6, double lambda_x0_8, double inv_gamma, double inv_gamma_2, double inv_gamma_3, double inv_gamma_4){
// //     double term1 = 6.5625 * inv_gamma_4;
// //     double term2 = 52.5 * lambda_x0_2 * inv_gamma_3;
// //     double term3 = 52.5 * lambda_x0_4 * inv_gamma_2;
// //     double term4 = 14.0 * lambda_x0_6 * inv_gamma;

// //     return term1 + term2 + term3 + term4 + lambda_x0_8;
// // }

// void P8_vec_opt3(double* lambda_x0_2, double* lambda_x0_4, double* lambda_x0_6, double* lambda_x0_8, double* inv_gamma, double* inv_gamma_2, double* inv_gamma_3, double* inv_gamma_4, double* result, int n) {
//     // Compute 6.5625 * inv_gamma_4[i] for each i, and store in result
//     vdMulC(n, inv_gamma_4, 6.5625, result);

//     // Compute 52.5 * lambda_x0_2[i] * inv_gamma_3[i] for each i, and store in a temporary array
//     double* temp1 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_2, inv_gamma_3, temp1); // temp1[i] = lambda_x0_2[i] * inv_gamma_3[i]
//     vdMulC(n, temp1, 52.5, temp1);             // temp1[i] = 52.5 * temp1[i]

//     // Compute 52.5 * lambda_x0_4[i] * inv_gamma_2[i] for each i, and store in another temporary array
//     double* temp2 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_4, inv_gamma_2, temp2); // temp2[i] = lambda_x0_4[i] * inv_gamma_2[i]
//     vdMulC(n, temp2, 52.5, temp2);             // temp2[i] = 52.5 * temp2[i]

//     // Compute 14.0 * lambda_x0_6[i] * inv_gamma[i] for each i, and store in yet another temporary array
//     double* temp3 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_6, inv_gamma, temp3);   // temp3[i] = lambda_x0_6[i] * inv_gamma[i]
//     vdMulC(n, temp3, 14.0, temp3);              // temp3[i] = 14.0 * temp3[i]

//     // Add temp1[i], temp2[i], and temp3[i] to result[i]
//     vdAdd(n, result, temp1, result);
//     vdAdd(n, result, temp2, result);
//     vdAdd(n, result, temp3, result);

//     // Add lambda_x0_8[i] to each result[i]
//     vdAdd(n, lambda_x0_8, result, result);

//     // Free the temporary arrays
//     mkl_free(temp1);
//     mkl_free(temp2);
//     mkl_free(temp3);
// }

// // inline double P9_opt2(double lambda_x0, double lambda_x0_3, double lambda_x0_5, double lambda_x0_7, double lambda_x0_9, double inv_gamma, double inv_gamma_2, double inv_gamma_3, double inv_gamma_4){
// //     double term1 = 59.0625 * lambda_x0 * inv_gamma_4;
// //     double term2 = 157.5 * lambda_x0_3 * inv_gamma_3;
// //     double term3 = 94.5 * lambda_x0_5 * inv_gamma_2;
// //     double term4 = 18.0 * lambda_x0_7 * inv_gamma;
// //     return term1 + term2 + term3 + term4 + lambda_x0_9;
// // }

// void P9_vec_opt3(double* lambda_x0, double* lambda_x0_3, double* lambda_x0_5, double* lambda_x0_7, double* lambda_x0_9, double* inv_gamma, double* inv_gamma_2, double* inv_gamma_3, double* inv_gamma_4, double* result, int n) {
//     // Compute 59.0625 * lambda_x0[i] * inv_gamma_4[i] for each i, and store in result
//     vdMul(n, lambda_x0, inv_gamma_4, result);   // result[i] = lambda_x0[i] * inv_gamma_4[i]
//     vdMulC(n, result, 59.0625, result);         // result[i] = 59.0625 * result[i]

//     // Compute 157.5 * lambda_x0_3[i] * inv_gamma_3[i] for each i, and store in a temporary array
//     double* temp1 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_3, inv_gamma_3, temp1);  // temp1[i] = lambda_x0_3[i] * inv_gamma_3[i]
//     vdMulC(n, temp1, 157.5, temp1);             // temp1[i] = 157.5 * temp1[i]

//     // Compute 94.5 * lambda_x0_5[i] * inv_gamma_2[i] for each i, and store in another temporary array
//     double* temp2 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_5, inv_gamma_2, temp2);  // temp2[i] = lambda_x0_5[i] * inv_gamma_2[i]
//     vdMulC(n, temp2, 94.5, temp2);              // temp2[i] = 94.5 * temp2[i]

//     // Compute 18.0 * lambda_x0_7[i] * inv_gamma[i] for each i, and store in yet another temporary array
//     double* temp3 = (double*)mkl_malloc(n * sizeof(double), 64); // 64-byte alignment for better performance with MKL
//     vdMul(n, lambda_x0_7, inv_gamma, temp3);    // temp3[i] = lambda_x0_7[i] * inv_gamma[i]
//     vdMulC(n, temp3, 18.0, temp3);               // temp3[i] = 18.0 * temp3[i]

//     // Add temp1[i], temp2[i], and temp3[i] to result[i]
//     vdAdd(n, result, temp1, result);
//     vdAdd(n, result, temp2, result);
//     vdAdd(n, result, temp3, result);

//     // Add lambda_x0_9[i] to each result[i]
//     vdAdd(n, lambda_x0_9, result, result);

//     // Free the temporary arrays
//     mkl_free(temp1);
//     mkl_free(temp2);
//     mkl_free(temp3);
// }


double calc_solid_MCSH_n1_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n)
{
    // value[0] +=  temp;
    double* temp = (double*)mkl_malloc(n * sizeof(double), 64);
    calc_temp_opt3(C1, C2, r0_sqr, occ, temp, n);

    double v_sum = 0;
    for (int i = 0; i < n; ++i) {
        v_sum += temp[i];
    }

    mkl_free(temp);
    return v_sum;
}

double calc_solid_MCSH_0_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n)
{

    // value[0] += temp;
    double* temp = (double*)mkl_malloc(n * sizeof(double), 64);
    calc_temp_opt3(C1, C2, r0_sqr, occ, temp, n);

    double v_sum = 0;
    for (int i = 0; i < n; ++i) {
        v_sum += temp[i];
    }

    mkl_free(temp);
    return v_sum;
}

double calc_solid_MCSH_1_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n)
{
    double* temp = (double*)mkl_malloc(n * sizeof(double), 64);
    calc_temp_opt3(C1, C2, r0_sqr, occ, temp, n);


    // double temp_x = lambda * x0;
    // double temp_y = lambda * y0;
    // double temp_z = lambda * z0;

    // double miu_1_1_1 = temp * temp_x;
    // double miu_1_1_2 = temp * temp_y;
    // double miu_1_1_3 = temp * temp_z;

    // value[0] += miu_1_1_1;
    // value[1] += miu_1_1_2;
    // value[2] += miu_1_1_3;

    double* temp_x = (double*)mkl_malloc(n * sizeof(double), 64);
    double* temp_y = (double*)mkl_malloc(n * sizeof(double), 64);
    double* temp_z = (double*)mkl_malloc(n * sizeof(double), 64);

    vdMul(n, lambda, x0, temp_x);
    vdMul(n, lambda, y0, temp_y);
    vdMul(n, lambda, z0, temp_z);

    vdMul(n, temp, temp_x, temp_x);
    vdMul(n, temp, temp_y, temp_y);
    vdMul(n, temp, temp_z, temp_z);

    double v0_sum = 0, v1_sum = 0, v2_sum = 0;
    for (int i = 0; i < n; ++i) {
        v0_sum += temp_x[i];
        v1_sum += temp_y[i];
        v2_sum += temp_z[i];
    }

    mkl_free(temp_x);
    mkl_free(temp_y);
    mkl_free(temp_z);

    mkl_free(temp);

    return v0_sum * v0_sum + v1_sum * v1_sum + v2_sum * v2_sum;

}


double calc_solid_MCSH_2_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n)
{
    double* temp = (double*)mkl_malloc(n * sizeof(double), 64);
    calc_temp_opt3(C1, C2, r0_sqr, occ, temp, n);

    // double lambda_x0 = lambda * x0;
    // double lambda_x0_2 = lambda_x0   * lambda_x0;
    // double lambda_y0 = lambda * y0;
    // double lambda_y0_2 = lambda_y0   * lambda_y0;
    // double lambda_z0 = lambda * z0;
    // double lambda_z0_2 = lambda_z0   * lambda_z0;
    // double inv_gamma = 1.0 / gamma;

    // double P2x = P2_opt2(lambda_x0_2, inv_gamma);
    // double P2y = P2_opt2(lambda_y0_2, inv_gamma);
    // double P2z = P2_opt2(lambda_z0_2, inv_gamma);

    double* lambda_x0 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_y0 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_z0 = (double*)mkl_malloc(n * sizeof(double), 64);

    vdMul(n, lambda, x0, lambda_x0);
    vdMul(n, lambda, y0, lambda_y0);
    vdMul(n, lambda, z0, lambda_z0);

    double* lambda_x0_2 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_y0_2 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_z0_2 = (double*)mkl_malloc(n * sizeof(double), 64);

    vdMul(n, lambda_x0, x0, lambda_x0_2);
    vdMul(n, lambda_y0, y0, lambda_y0_2);
    vdMul(n, lambda_z0, z0, lambda_z0_2);

    double* inv_gamma = (double*)mkl_malloc(n * sizeof(double), 64);

    vdInv(n, gamma, inv_gamma);

    double* P2x = (double*)mkl_malloc(n * sizeof(double), 64);
    double* P2y = (double*)mkl_malloc(n * sizeof(double), 64);
    double* P2z = (double*)mkl_malloc(n * sizeof(double), 64);

    P2_vec_opt3(lambda_x0_2, inv_gamma, P2x, n);
    P2_vec_opt3(lambda_y0_2, inv_gamma, P2y, n);
    P2_vec_opt3(lambda_z0_2, inv_gamma, P2z, n);

    // group 1

    // double gp1_term_x = (2.0 * P2x) - (P2y + P2z);
    // double gp1_term_y = (2.0 * P2y) - (P2x + P2z);
    // double gp1_term_z = (2.0 * P2z) - (P2x + P2y);

    // double gp1_miu_1 = temp * gp1_term_x;
    // double gp1_miu_2 = temp * gp1_term_y;
    // double gp1_miu_3 = temp * gp1_term_z;

    // value[0] += gp1_miu_1;
    // value[1] += gp1_miu_2;
    // value[2] += gp1_miu_3;

    double* v0 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v1 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v2 = (double*)mkl_malloc(n * sizeof(double), 64);

    vdAdd(n, P2y, P2z, v0);
    cblas_daxpby(n, 2.0, P2x, 1, -1.0, v0, 1);
    vdMul(n, v0, temp, v0); 

    vdAdd(n, P2x, P2z, v1);
    cblas_daxpby(n, 2.0, P2y, 1, -1.0, v1, 1);
    vdMul(n, v1, temp, v1); 

    vdAdd(n, P2x, P2y, v2);
    cblas_daxpby(n, 2.0, P2z, 1, -1.0, v2, 1);
    vdMul(n, v2, temp, v2); 



    // group 2

    // double gp2_miu_2_2_1 = 3.0 * temp * lambda_x0 * lambda_y0;
    // double gp2_miu_2_2_2 = 3.0 * temp * lambda_x0 * lambda_z0;
    // double gp2_miu_2_2_3 = 3.0 * temp * lambda_y0 * lambda_z0;

    // value[3] += gp2_miu_2_2_1;
    // value[4] += gp2_miu_2_2_2;
    // value[5] += gp2_miu_2_2_3;

    double* v3 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v4 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v5 = (double*)mkl_malloc(n * sizeof(double), 64);

    vdMul(n, lambda_x0, lambda_y0, v3);
    vdMul(n, temp, v3, v3);
    cblas_dscal(n, 3.0, v3, 1);

    vdMul(n, lambda_x0, lambda_z0, v4);
    vdMul(n, temp, v4, v4);
    cblas_dscal(n, 3.0, v4, 1);

    vdMul(n, lambda_y0, lambda_z0, v5);
    vdMul(n, temp, v5, v5);
    cblas_dscal(n, 3.0, v5, 1);



    double v0_sum = 0, v1_sum = 0, v2_sum = 0, v3_sum = 0, v4_sum = 0, v5_sum = 0;
    for (int i = 0; i < n; ++i) {
        v0_sum += v0[i];
        v1_sum += v1[i];
        v2_sum += v2[i];
        v3_sum += v3[i];
        v4_sum += v4[i];
        v5_sum += v5[i];
    }

    mkl_free(lambda_x0);
    mkl_free(lambda_y0);
    mkl_free(lambda_z0);
    mkl_free(lambda_x0_2);
    mkl_free(lambda_y0_2);
    mkl_free(lambda_z0_2);
    mkl_free(inv_gamma);
    mkl_free(P2x);
    mkl_free(P2y);
    mkl_free(P2z);

    mkl_free(v0);
    mkl_free(v1);
    mkl_free(v2);
    mkl_free(v3);
    mkl_free(v4);
    mkl_free(v5);

    mkl_free(temp);

    return (v0_sum * v0_sum + v1_sum * v1_sum + v2_sum * v2_sum) + 2.0 * (v3_sum * v3_sum + v4_sum * v4_sum + v5_sum * v5_sum);

}

double calc_solid_MCSH_3_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n)
{
    double* temp = (double*)mkl_malloc(n * sizeof(double), 64);
    calc_temp_opt3(C1, C2, r0_sqr, occ, temp, n);

    // double lambda_x0 = lambda * x0;
    // double lambda_x0_2 = lambda_x0   * lambda_x0;
    // double lambda_x0_3 = lambda_x0_2 * lambda_x0;
    // double lambda_y0 = lambda * y0;
    // double lambda_y0_2 = lambda_y0   * lambda_y0;
    // double lambda_y0_3 = lambda_y0_2 * lambda_y0;
    // double lambda_z0 = lambda * z0;
    // double lambda_z0_2 = lambda_z0   * lambda_z0;
    // double lambda_z0_3 = lambda_z0_2 * lambda_z0;

    // double inv_gamma = 1.0 / gamma;

    // double P1x = lambda_x0;
    // double P1y = lambda_y0;
    // double P1z = lambda_z0;

    // double P2x = P2_opt2(lambda_x0_2, inv_gamma);
    // double P2y = P2_opt2(lambda_y0_2, inv_gamma);
    // double P2z = P2_opt2(lambda_z0_2, inv_gamma);

    // double P3x = P3_opt2(lambda_x0, lambda_x0_3, inv_gamma);
    // double P3y = P3_opt2(lambda_y0, lambda_y0_3, inv_gamma);
    // double P3z = P3_opt2(lambda_z0, lambda_z0_3, inv_gamma);


    double* lambda_x0 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_y0 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_z0 = (double*)mkl_malloc(n * sizeof(double), 64);

    vdMul(n, lambda, x0, lambda_x0);
    vdMul(n, lambda, y0, lambda_y0);
    vdMul(n, lambda, z0, lambda_z0);

    double* lambda_x0_2 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_y0_2 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_z0_2 = (double*)mkl_malloc(n * sizeof(double), 64);

    vdMul(n, lambda_x0, x0, lambda_x0_2);
    vdMul(n, lambda_y0, y0, lambda_y0_2);
    vdMul(n, lambda_z0, z0, lambda_z0_2);

    double* lambda_x0_3 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_y0_3 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* lambda_z0_3 = (double*)mkl_malloc(n * sizeof(double), 64);

    vdMul(n, lambda_x0_2, x0, lambda_x0_3);
    vdMul(n, lambda_y0_2, y0, lambda_y0_3);
    vdMul(n, lambda_z0_2, z0, lambda_z0_3);

    double* inv_gamma = (double*)mkl_malloc(n * sizeof(double), 64);

    vdInv(n, gamma, inv_gamma);

    double* P2x = (double*)mkl_malloc(n * sizeof(double), 64);
    double* P2y = (double*)mkl_malloc(n * sizeof(double), 64);
    double* P2z = (double*)mkl_malloc(n * sizeof(double), 64);

    P2_vec_opt3(lambda_x0_2, inv_gamma, P2x, n);
    P2_vec_opt3(lambda_y0_2, inv_gamma, P2y, n);
    P2_vec_opt3(lambda_z0_2, inv_gamma, P2z, n);

    double* P3x = (double*)mkl_malloc(n * sizeof(double), 64);
    double* P3y = (double*)mkl_malloc(n * sizeof(double), 64);
    double* P3z = (double*)mkl_malloc(n * sizeof(double), 64);

    P3_vec_opt3(lambda_x0, lambda_x0_3, inv_gamma, P3x, n);
    P3_vec_opt3(lambda_x0, lambda_y0_3, inv_gamma, P3y, n);
    P3_vec_opt3(lambda_x0, lambda_z0_3, inv_gamma, P3z, n);

    // group 1

    // double gp1_term_x = (6.0 * P3x) - (9.0 * P1x * (P2y + P2z));
    // double gp1_term_y = (6.0 * P3y) - (9.0 * P1y * (P2x + P2z));
    // double gp1_term_z = (6.0 * P3z) - (9.0 * P1z * (P2x + P2y));

    // double gp1_miu_1 = temp * gp1_term_x;
    // double gp1_miu_2 = temp * gp1_term_y;
    // double gp1_miu_3 = temp * gp1_term_z;

    // value[0] += gp1_miu_1;
    // value[1] += gp1_miu_2;
    // value[2] += gp1_miu_3;


    double* v0 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v1 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v2 = (double*)mkl_malloc(n * sizeof(double), 64);


    vdAdd(n, P2y, P2z, v0);
    vdMul(n, lambda_x0, v0, v0);
    cblas_daxpby(n, 6.0, P3x, 1, -9.0, v0, 1);
    vdMul(n, v0, temp, v0);

    vdAdd(n, P2x, P2z, v1);
    vdMul(n, lambda_y0, v1, v1);
    cblas_daxpby(n, 6.0, P3y, 1, -9.0, v1, 1);
    vdMul(n, v1, temp, v1);

    vdAdd(n, P2x, P2y, v2);
    vdMul(n, lambda_z0, v2, v2);
    cblas_daxpby(n, 6.0, P3z, 1, -9.0, v2, 1);
    vdMul(n, v2, temp, v2);


    // group 2

    // double gp2_term_1 = (12.0 * P1y * P2x) - (3.0 * P3y) - (3.0 * P1y * P2z);
    // double gp2_term_2 = (12.0 * P1x * P2y) - (3.0 * P3x) - (3.0 * P1x * P2z);
    // double gp2_term_3 = (12.0 * P1z * P2x) - (3.0 * P3z) - (3.0 * P1z * P2y);
    // double gp2_term_4 = (12.0 * P1x * P2z) - (3.0 * P3x) - (3.0 * P1x * P2y);
    // double gp2_term_5 = (12.0 * P1z * P2y) - (3.0 * P3z) - (3.0 * P1z * P2x);
    // double gp2_term_6 = (12.0 * P1y * P2z) - (3.0 * P3y) - (3.0 * P1y * P2x);

    // double gp2_miu_1 = temp * gp2_term_1;
    // double gp2_miu_2 = temp * gp2_term_2;
    // double gp2_miu_3 = temp * gp2_term_3;
    // double gp2_miu_4 = temp * gp2_term_4;
    // double gp2_miu_5 = temp * gp2_term_5;
    // double gp2_miu_6 = temp * gp2_term_6;

    // value[3] += gp2_miu_1;
    // value[4] += gp2_miu_2;
    // value[5] += gp2_miu_3;
    // value[6] += gp2_miu_4;
    // value[7] += gp2_miu_5;
    // value[8] += gp2_miu_6;


    double* gp2_temp_xy = (double*)mkl_malloc(n * sizeof(double), 64);
    double* gp2_temp_xz = (double*)mkl_malloc(n * sizeof(double), 64);
    double* gp2_temp_yz = (double*)mkl_malloc(n * sizeof(double), 64);
    double* gp2_temp_yx = (double*)mkl_malloc(n * sizeof(double), 64);
    double* gp2_temp_zx = (double*)mkl_malloc(n * sizeof(double), 64);
    double* gp2_temp_zy = (double*)mkl_malloc(n * sizeof(double), 64);

    vdMul(n, lambda_x0, P2y, gp2_temp_xy);
    vdMul(n, lambda_x0, P2z, gp2_temp_xz);
    vdMul(n, lambda_y0, P2z, gp2_temp_yz);
    vdMul(n, lambda_y0, P2x, gp2_temp_yx);
    vdMul(n, lambda_z0, P2x, gp2_temp_zx);
    vdMul(n, lambda_z0, P2y, gp2_temp_zy);

    double* v3 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v4 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v5 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v6 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v7 = (double*)mkl_malloc(n * sizeof(double), 64);
    double* v8 = (double*)mkl_malloc(n * sizeof(double), 64);

    cblas_dcopy(n, gp2_temp_yz, 1, v3, 1);
    cblas_dcopy(n, gp2_temp_xz, 1, v4, 1);
    cblas_dcopy(n, gp2_temp_zy, 1, v5, 1);
    cblas_dcopy(n, gp2_temp_xy, 1, v6, 1);
    cblas_dcopy(n, gp2_temp_zx, 1, v7, 1);
    cblas_dcopy(n, gp2_temp_yx, 1, v8, 1);

    cblas_daxpby(n, 12.0, gp2_temp_yx, 1, -3.0, v3, 1);
    cblas_daxpby(n, 12.0, gp2_temp_xy, 1, -3.0, v4, 1);
    cblas_daxpby(n, 12.0, gp2_temp_zx, 1, -3.0, v5, 1);
    cblas_daxpby(n, 12.0, gp2_temp_xz, 1, -3.0, v6, 1);
    cblas_daxpby(n, 12.0, gp2_temp_zy, 1, -3.0, v7, 1);
    cblas_daxpby(n, 12.0, gp2_temp_yz, 1, -3.0, v8, 1);

    cblas_daxpy(n, -3.0, P3y, 1, v3, 1);
    cblas_daxpy(n, -3.0, P3x, 1, v4, 1);
    cblas_daxpy(n, -3.0, P3z, 1, v5, 1);
    cblas_daxpy(n, -3.0, P3x, 1, v6, 1);
    cblas_daxpy(n, -3.0, P3z, 1, v7, 1);
    cblas_daxpy(n, -3.0, P3y, 1, v8, 1);


    vdMul(n, v3, temp, v3);
    vdMul(n, v4, temp, v4);
    vdMul(n, v5, temp, v5);
    vdMul(n, v6, temp, v6);
    vdMul(n, v7, temp, v7);
    vdMul(n, v8, temp, v8);

    mkl_free(gp2_temp_xy);
    mkl_free(gp2_temp_yx);
    mkl_free(gp2_temp_xz);
    mkl_free(gp2_temp_zx);
    mkl_free(gp2_temp_yz);
    mkl_free(gp2_temp_zy);

    // group 3

    // double gp3_m_3_3 = 15.0 * temp * lambda_x0 * lambda_y0 * lambda_z0; 

    // value[9] += gp3_m_3_3;

    double* v9 = (double*)mkl_malloc(n * sizeof(double), 64);

    // vdMulC(n, temp, 15.0, v9);
    vdMul(n, temp, lambda_x0, v9);
    vdMul(n, v9, lambda_y0, v9);
    vdMul(n, v9, lambda_z0, v9);
    cblas_dscal(n, 15.0, v9, 1);


    double v0_sum = 0, v1_sum = 0, v2_sum = 0, v3_sum = 0, v4_sum = 0, v5_sum = 0, v6_sum = 0, v7_sum = 0, v8_sum = 0, v9_sum = 0;
    for (int i = 0; i < n; ++i) {
        v0_sum += v0[i];
        v1_sum += v1[i];
        v2_sum += v2[i];
        v3_sum += v3[i];
        v4_sum += v4[i];
        v5_sum += v5[i];
        v6_sum += v6[i];
        v7_sum += v7[i];
        v8_sum += v8[i];
        v9_sum += v9[i];
    }


    mkl_free(lambda_x0);
    mkl_free(lambda_y0);
    mkl_free(lambda_z0);
    mkl_free(lambda_x0_2);
    mkl_free(lambda_y0_2);
    mkl_free(lambda_z0_2);
    mkl_free(lambda_x0_3);
    mkl_free(lambda_y0_3);
    mkl_free(lambda_z0_3);
    mkl_free(inv_gamma);
    mkl_free(P2x);
    mkl_free(P2y);
    mkl_free(P2z);
    mkl_free(P3x);
    mkl_free(P3y);
    mkl_free(P3z);

    mkl_free(v0);
    mkl_free(v1);
    mkl_free(v2);
    mkl_free(v3);
    mkl_free(v4);
    mkl_free(v5);
    mkl_free(v6);
    mkl_free(v7);
    mkl_free(v8);
    mkl_free(v9);

    mkl_free(temp);

    return (v0_sum * v0_sum + v1_sum * v1_sum + v2_sum * v2_sum) 
            + 3.0 * (v3_sum * v3_sum + v4_sum * v4_sum + v5_sum * v5_sum + v6_sum * v6_sum + v7_sum * v7_sum + v8_sum * v8_sum)
            + 6.0 * v9_sum * v9_sum;
}












