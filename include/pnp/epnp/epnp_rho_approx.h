
#pragma once
#define MAT_API extern
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "epnp_data.h"
#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus


MAT_API void  copy_R_and_t(const double R_src[3][3], const double t_src[3],
                          double R_dst[3][3], double t_dst[3]);
MAT_API double  dist2(const double * p1, const double * p2);
MAT_API double  dot(const double * v1, const double * v2);
MAT_API void  estimate_R_and_t(double R[3][3], double t[3],epnp_data *d);
MAT_API void compute_ccs(const double * betas, const double * ut,epnp_data *d);
MAT_API void compute_pcs(epnp_data *d);
MAT_API void  solve_for_sign(epnp_data *d);
MAT_API double  compute_R_and_t(const double * ut, const double * betas,
                               double R[3][3], double t[3],epnp_data *d);
MAT_API void  find_betas_approx_1(const double * L_6x10, const double * Rho,
                                 double * betas);
MAT_API void  find_betas_approx_2(const double * L_6x10, const double * Rho,
                                 double * betas);
MAT_API void  find_betas_approx_3(const double * L_6x10, const double * Rho,
                                 double * betas);
MAT_API void  compute_L_6x10(const double * ut, double * l_6x10);
MAT_API void  compute_rho(double * rho ,epnp_data *d);
MAT_API void  compute_A_and_b_gapoint_2d_arrays_newton(const double * l_6x10, const double * rho,
                                                      double betas[4], double * A, double * b);
MAT_API void  qr_solve(double * A, double * b, double * X,int m,int n);
MAT_API void  gapoint_2d_arrays_newton(const double * L_6x10, const double * Rho,
                                      double betas[4]);

MAT_API double   reprojection_error(epnp_data *d, double R[3][3],  double t[3]);

#ifdef __cplusplus
}
#endif // __cplusplus
