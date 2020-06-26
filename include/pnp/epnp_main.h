
#pragma once
#define MAT_API extern
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "epnp/epnp_data.h"
typedef struct Point_double
{
    double x;
    double y;
}Point_double;

typedef struct Point3d_double
{
    double x;
    double y;
    double z;
}Point3d_double;

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

MAT_API double epnp_computePose_main(epnp_data * PnP ,float *answer,int width,int height,int fovX,int fovY,Point_double* finded, Point3d_double *knownPos);
MAT_API void init_test_points(Point_double* finded);
MAT_API void epnp_init(epnp_data * PnP,int width,int height,int fovX,int fovY);

inline MAT_API void buildProjectionMatrix(const epnp_data *c, float *out_mat ,float far,float near);


#ifdef __cplusplus
}
#endif // __cplusplus



