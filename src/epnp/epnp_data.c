
#include <pnp/epnp/epnp_data.h>
#include <pnp/epnp/math/matrix.h>

#include "math.h"

MAT_API void set_maximum_number_of_correspondences(int n,epnp_data* d)
{

    d->maximum_number_of_correspondences = n;
}




MAT_API void  set_internal_parameters(epnp_data * d,int width, int height, int fovX, int fovY)
{
    d->width=width;
    d->height=height;
    d->center_x = (double)width/2.0;
    d->center_y = (double)height/2.0;
    d->focal_lenght_x = d->center_x/tan( (   (double)fovX/180.)*M_PI );

    d->focal_lenght_y = d->center_y/tan( (   (double)fovY/180.)*M_PI );
}

MAT_API void  reset_correspondences(epnp_data* d)
{

    d->number_of_correspondences = 0;
}

MAT_API void  add_correspondence(epnp_data *d,double X, double Y, double Z, double u, double v)
{
    d->point_3d_array[3 * d->number_of_correspondences    ] = X;
    d->point_3d_array[3 * d->number_of_correspondences + 1] = Y;
    d->point_3d_array[3 * d->number_of_correspondences + 2] = Z;


    d->point_2d_array[2 * d->number_of_correspondences    ] = u;
    d->point_2d_array[2 * d->number_of_correspondences + 1] = v;

    d->number_of_correspondences++;
}




//start
MAT_API void prepare_data(epnp_data* d)
{
    // d->cws[0][0] = d->cws[0][1] = d->cws[0][2] = 0;

    int var=0;
    int var2=0;
    for (var = 0; var < 4; ++var)
    {
        for (var2 = 0; var2 < 3; ++var2)
        {
            d->cws[var][var2]=0.0;
        }
    }

    int i=0;
    int j=0;

    for( i = 0; i < d->number_of_correspondences; i++)
    {
        for( j = 0; j < 3; j++)
        {
            d->cws[0][j] += d->point_3d_array[3 * i + j];

        }
    }



    for(int j = 0; j < 3; j++)
        d->cws[0][j] /= d->number_of_correspondences;



}



