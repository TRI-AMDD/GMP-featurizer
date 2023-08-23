#include "mkl.h"

typedef double (*SolidGMPFunctionNoderivOpt3) (double*, double*, double*, double*, double*, double*, double*, double*, double*, int);
SolidGMPFunctionNoderivOpt3 get_solid_mcsh_function_noderiv_opt3(int mcsh_order);

void calc_temp_opt3(double* C1, double* C2, double* r0_sqr, double* occ, double* result, int n);

void P2_vec_opt3(double* lambda_x0_2, double* inv_gamma, double* result, int n);
void P3_vec_opt3(double* lambda_x0, double* lambda_x0_3, double* inv_gamma, double* result, int n);
void P4_vec_opt3(double* lambda_x0_2, double* lambda_x0_4, double* inv_gamma, double* inv_gamma_2, double* result, int n);
void P5_vec_opt3(double* lambda_x0, double* lambda_x0_3, double* lambda_x0_5, double* inv_gamma, double* inv_gamma_2, double* result, int n);
void P6_vec_opt3(double* lambda_x0_2, double* lambda_x0_4, double* lambda_x0_6, double* inv_gamma, double* inv_gamma_2, double* inv_gamma_3, double* result, int n);
void P7_vec_opt3(double* lambda_x0, double* lambda_x0_3, double* lambda_x0_5, double* lambda_x0_7, double* inv_gamma, double* inv_gamma_2, double* inv_gamma_3, double* result, int n);
void P8_vec_opt3(double* lambda_x0_2, double* lambda_x0_4, double* lambda_x0_6, double* lambda_x0_8, double* inv_gamma, double* inv_gamma_2, double* inv_gamma_3, double* inv_gamma_4, double* result, int n);
void P9_vec_opt3(double* lambda_x0, double* lambda_x0_3, double* lambda_x0_5, double* lambda_x0_7, double* lambda_x0_9, double* inv_gamma, double* inv_gamma_2, double* inv_gamma_3, double* inv_gamma_4, double* result, int n);

double calc_solid_MCSH_n1_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n);
double calc_solid_MCSH_0_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n);
double calc_solid_MCSH_1_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n);
double calc_solid_MCSH_2_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n);
double calc_solid_MCSH_3_noderiv_opt3(double* x0, double* y0, double* z0, double* r0_sqr, double* occ, double* C1, double* C2, double* lambda, double* gamma, int n);
// double calc_solid_MCSH_4_noderiv_opt3(double* x0, double* y0, double* z0, double* temp, double* lambda, double* gamma, int n);
