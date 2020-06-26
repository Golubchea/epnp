
#pragma once
#define MAT_API extern
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
typedef struct epnp_data
{
    //in params
    int width;
    int height;

    double center_x;
    double center_y;

    double  focal_lenght_x;
    double  focal_lenght_y;
    unsigned int maximum_number_of_correspondences;
    //

    unsigned int n;

    int number_of_correspondences;
    double cws[4][3], ccs[4][3];

    double point_3d_array[3*4];
    double point_2d_array[2*4];
    double alphas_array[4*4];
    double pcs_array[3*4];

}epnp_data;

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

MAT_API void  set_internal_parameters(epnp_data * d,int width, int height, int fovX, int fovY);
MAT_API void  set_maximum_number_of_correspondences(int n,epnp_data* d);
MAT_API void  reset_correspondences(epnp_data* d);
MAT_API void  add_correspondence(epnp_data *d,double X, double Y, double Z, double u, double v);

MAT_API void prepare_data(epnp_data* d);


#ifdef __cplusplus
}
#endif // __cplusplus





