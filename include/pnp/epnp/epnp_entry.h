
#pragma once
#define MAT_API extern
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "epnp_data.h"
#include "math/singular_value_decomposition.h"
#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

MAT_API double  compute_pose(epnp_data *PnP,double R_est [3][3],double t_est[3]);
MAT_API void  compute_barycentric_coordinates(epnp_data* d);

MAT_API void  getModelMat(double R[3][3], double t[3],float *out);
#ifdef __cplusplus
}
#endif // __cplusplus
double  compute_pose2(epnp_data *PnP,double R_est [3][3],double t_est[3]);
