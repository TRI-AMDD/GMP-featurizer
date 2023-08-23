
#include "helper_optimization2.h"

// double P1_opt2(double lambda_x0){
//     return lambda_x0;
// }

// double P2_opt2(double lambda_x0_2, double inv_gamma){
//     // double lambda_x0_2 = lambda * lambda * x0 * x0;
//     return (0.5 * inv_gamma) + lambda_x0_2;
// }

// double P3_opt2(double lambda_x0, double lambda_x0_3, double inv_gamma){
//     // double lambda_x0 = lambda * x0;
//     // double lambda_x0_3 = lambda_x0 * lambda_x0 * lambda_x0;
//     return (1.5 * lambda_x0 * inv_gamma) + lambda_x0_3;
// }

// double P4_opt2(double lambda_x0_2, double lambda_x0_4, double inv_gamma, double inv_gamma_2){
//     // double lambda_x0_2 = lambda * lambda * x0 * x0;
//     // double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
//     return (0.75 * inv_gamma_2) + (3.0 * lambda_x0_2 * inv_gamma) + lambda_x0_4;
// }

// double P5_opt2(double lambda_x0, double lambda_x0_3, double lambda_x0_5, double inv_gamma, double inv_gamma_2){
//     // double lambda_x0 = lambda * x0;
//     // double lambda_x0_2 = lambda_x0 * lambda_x0;
//     // double lambda_x0_3 = lambda_x0 * lambda_x0_2;
//     // double lambda_x0_5 = lambda_x0_3 * lambda_x0_2;
//     return ((15.0 * lambda_x0) * (0.25 * inv_gamma_2)) + (5.0 * lambda_x0_3 * inv_gamma) + lambda_x0_5;
// }

// double P6_opt2(double lambda_x0_2, double lambda_x0_4, double lambda_x0_6, double inv_gamma, double inv_gamma_2, double inv_gamma_3){
//     // double lambda_x0_2 = lambda * lambda * x0 * x0;
//     // double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
//     // double lambda_x0_6 = lambda_x0_4 * lambda_x0_2;
//     return (1.875 * inv_gamma_3) + (11.25 * lambda_x0_2 * inv_gamma_2) + (7.5 * lambda_x0_4 *inv_gamma) + lambda_x0_6;
// }

// double P7_opt2(double lambda_x0, double lambda_x0_3, double lambda_x0_5, double lambda_x0_7, double inv_gamma, double inv_gamma_2, double inv_gamma_3){
//     // double lambda_x0 = lambda * x0;
//     // double lambda_x0_2 = lambda_x0 * lambda_x0;
//     // double lambda_x0_3 = lambda_x0 * lambda_x0_2;
//     // double lambda_x0_5 = lambda_x0_3 * lambda_x0_2;
//     // double lambda_x0_7 = lambda_x0_5 * lambda_x0_2;
//     double term1 = 13.125 * lambda_x0 * inv_gamma_3;
//     double term2 = 26.25 * lambda_x0_3 * inv_gamma_2;
//     double term3 = 10.5 * lambda_x0_5 * inv_gamma;
//     return term1 + term2 + term3 + lambda_x0_7;
// }

// double P8_opt2(double lambda_x0_2, double lambda_x0_4, double lambda_x0_6, double lambda_x0_8, double inv_gamma, double inv_gamma_2, double inv_gamma_3, double inv_gamma_4){
//     // double lambda_x0_2 = lambda * lambda * x0 * x0;
//     // double lambda_x0_4 = lambda_x0_2 * lambda_x0_2;
//     // double lambda_x0_6 = lambda_x0_4 * lambda_x0_2;
//     // double lambda_x0_8 = lambda_x0_6 * lambda_x0_2;
//     double term1 = 6.5625 * inv_gamma_4;
//     double term2 = 52.5 * lambda_x0_2 * inv_gamma_3;
//     double term3 = 52.5 * lambda_x0_4 * inv_gamma_2;
//     double term4 = 14.0 * lambda_x0_6 * inv_gamma;

//     return term1 + term2 + term3 + term4 + lambda_x0_8;
// }

// double P9_opt2(double lambda_x0, double lambda_x0_3, double lambda_x0_5, double lambda_x0_7, double lambda_x0_9, double inv_gamma, double inv_gamma_2, double inv_gamma_3, double inv_gamma_4){
//     // double lambda_x0 = lambda * x0;
//     // double lambda_x0_2 = lambda_x0 * lambda_x0;
//     // double lambda_x0_3 = lambda_x0 * lambda_x0_2;
//     // double lambda_x0_5 = lambda_x0_3 * lambda_x0_2;
//     // double lambda_x0_7 = lambda_x0_5 * lambda_x0_2;
//     // double lambda_x0_9 = lambda_x0_7 * lambda_x0_2;
//     double term1 = 59.0625 * lambda_x0 * inv_gamma_4;
//     double term2 = 157.5 * lambda_x0_3 * inv_gamma_3;
//     double term3 = 94.5 * lambda_x0_5 * inv_gamma_2;
//     double term4 = 18.0 * lambda_x0_7 * inv_gamma;
//     // double term1 = 945.0 * lambda_x0 / (16.0 * gamma * gamma * gamma * gamma);
//     // double term2 = 315.0 * lambda_x0_3 / (2.0 * gamma * gamma * gamma);
//     // double term3 = 189.0 * lambda_x0_5 / (2.0 * gamma * gamma);
//     // double term4 = 18.0 * lambda_x0_7 / gamma;
//     return term1 + term2 + term3 + term4 + lambda_x0_9;
// }


