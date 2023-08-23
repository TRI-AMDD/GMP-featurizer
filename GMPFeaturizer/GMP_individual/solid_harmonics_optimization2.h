#include "helper.h"

typedef void (*SolidGMPFunctionNoderivOpt2) (double, double, double, double, double, double, double, double *);
SolidGMPFunctionNoderivOpt2 get_solid_mcsh_function_noderiv_opt2(int mcsh_order);

void calc_solid_MCSH_n1_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_0_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_1_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_2_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_3_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_4_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_5_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_6_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_7_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_8_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);
void calc_solid_MCSH_9_noderiv_opt2(double x0, double y0, double z0, double r0_sqr, double temp, double lambda, double gamma, double *value);

double get_desc_value_opt2(int mcsh_order, double* v);
int get_num_order_values(int mcsh_order);