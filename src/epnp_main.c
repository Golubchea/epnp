#include <pnp/epnp_main.h>
#include <pnp/epnp/epnp_data.h>
#include <pnp/epnp/epnp_data.h>
#include <pnp/epnp/epnp_entry.h>
#include <pnp/containers/vector.h>


inline MAT_API void buildProjectionMatrix(const epnp_data *c, float *out_mat ,float far,float near)
{

    assert(c->focal_lenght_x > 0 );

    float s=0;


    out_mat[0]=2.0*c->focal_lenght_x/c->width;
    out_mat[1]= 0.;
    out_mat[2]= 0.;
    out_mat[3]= 0.;

    float t1=(far-near);

    out_mat[4]= 2.0f*s/(float)c->width;
    out_mat[5]= 2.0f*(float)c->focal_lenght_y/(float)c->height;
    out_mat[6]= 0.0 ;
    out_mat[7]= 0.0 ;

    out_mat[8]=  2.0f*((float)c->center_x/(float)c->width)-1.0f;
    out_mat[9]=  2.0f*((float)c->center_y/(float)c->height)-1.0f;
    out_mat[10]= -(far+near)/t1;
    out_mat[11]=  -1.0;

    out_mat[12]= 0;
    out_mat[13]= 0;
    out_mat[14]= -2.0f*far*near/t1;
    out_mat[15]= 0;

}

MAT_API void epnp_init(epnp_data * PnP,int width,int height,int fovX,int fovY)
{
    reset_correspondences(PnP);

    set_internal_parameters(PnP,width,height,fovX,fovY);


    set_maximum_number_of_correspondences(4,PnP);
}



MAT_API double epnp_computePose_main(epnp_data * PnP, float *answer,int width,int height,int fovX,int fovY,Point_double* finded,Point3d_double *known)
{

    reset_correspondences(PnP);
    if(width!=PnP->width ||  height!=PnP->height)
    {
        set_internal_parameters(PnP,width,height,fovX,fovY);
    }
    else
    {

    }

    set_maximum_number_of_correspondences(4,PnP);

    //4--2
    //|  |
    //3--1

    //for test
//    add_correspondence(PnP,1.0, 0.0, 0.0, finded[0].x, finded[0].y);
//    add_correspondence(PnP,1.0, 0.0, 1.0, finded[1].x, finded[1].y);
//    add_correspondence(PnP,0.0, 0.0, 0.0, finded[2].x, finded[2].y);
//    add_correspondence(PnP,0.0, 0.0, 1.0, finded[3].x, finded[3].y);
                              //x              z              y
    add_correspondence(PnP,known[0].x,        known[0].y,  known[0].z,   finded[0].x, finded[0].y);
    add_correspondence(PnP,known[1].x,        known[1].y,  known[1].z,   finded[1].x, finded[1].y);
    add_correspondence(PnP,known[2].x,        known[2].y,  known[2].z,   finded[2].x, finded[2].y);
    add_correspondence(PnP,known[3].x,        known[3].y,  known[3].z,   finded[3].x, finded[3].y);

    double R_est[3][3], t_est[3];

#ifdef __cplusplus
    double err2 =  compute_pose2(PnP,R_est, t_est);
#endif

#ifndef __cplusplus
    double Errortica =  compute_pose(PnP,R_est, t_est);
#endif




    //debug SIGN SHIT!!!!!!
    //caused by SVD decomposition is not unique
//    R_est[0][1]=-R_est[0][1];
//    R_est[1][1]=-R_est[1][1];
//    R_est[2][0]=-R_est[2][0];
//    R_est[2][2]=-R_est[2][2];

    getModelMat(R_est,t_est,answer);

    return Errortica;
}



MAT_API void init_test_points(Point_double* finded)
{
    finded[0].x=335.000000;  finded[0].y=321.00000;
    finded[1].x=448.000000;  finded[1].y=155.00000;
    finded[2].x=128.000000;  finded[2].y=207.00000;
    finded[3].x=263.000000;  finded[3].y=75.000000;


}
