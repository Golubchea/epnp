#include <pnp/epnp/epnp_data.h>
#include "math.h"

#include <pnp/epnp/math/matrix.h>
#include <pnp/epnp/math/singular_value_decomposition.h>
//part2

#include <pnp/epnp/epnp_rho_approx.h>
#include <stdlib.h>
#define getName(var)  #var

//#include "math/SVD.h"

void choose_control_points(epnp_data* d)
{
    prepare_data(d);


    double PW0[4*3];
    // Take C1, C2, and C3 from PCA on the reference points:

    for(int i = 0; i < d->number_of_correspondences; i++)
        for(int j = 0; j < 3; j++)
            PW0[3 * i + j] = d->point_3d_array[3 * i + j] - d->cws[0][j];

    double res[3*4];
    Atranspose_mul_A(PW0,res,4,3);



    //print_mMxN_1(res,3,3,getName(res));



    double U[3][3];
    double V[3][3];
    double singular_values[3];


    int error;
    error = Singular_Value_Decomposition(res, 3, 3, (double*) U,
                                         singular_values, (double*) V);


    matrix_transpose2( (double*)U,(double*)V,3,3);
    //  matrix_mul_scalar((double*)V,-1.0,3,3);


    // matrix_print(singular_values,1,3,getName(singular_values));
    // matrix_print((double*)V,3,3,getName(V));

    double *Uptr=(double*)V;
    for(int i = 1; i < 4; i++)
    {
        double k = sqrt( singular_values[i - 1] / d->number_of_correspondences);
        for(int j = 0; j < 3; j++)
        {
            d->cws[i][j] = d->cws[0][j] + k *  Uptr[3 * (i - 1) + j];

        }

    }
    // matrix_print((double*) d->cws,4,3,getName(d->cws));

}

MAT_API void  compute_barycentric_coordinates(epnp_data* d)
{
    double cc[3 * 3], cc_inv[3 * 3];

    for(int i = 0; i < 3; i++)
    {
        for(int j = 1; j < 4; j++)
        {
            cc[3 * i + j - 1] = d->cws[j][i] - d->cws[0][i];
        }
    }

    double U[3][3];
    double V[3][3];
    double singular_values[3];



    int error;
    error = Singular_Value_Decomposition(cc, 3, 3, (double*) V,
                                         singular_values, (double*) U);




    double Astar[3][3];
    double tolerance=0.00001;

    Singular_Value_Decomposition_Inverse((double*) U, singular_values, (double*) V,
                                         tolerance, 3, 3, (double*) Astar);



    matrix_transpose2(  (double*)Astar,(double*)cc_inv,3,3);


    // matrix_print(cc,3,3,"cc");

    // matrix_print((double*)cc_inv,3,3,"cc_inv");
    //    cvInvert(&CC, &CC_inv, CV_SVD);
    double * ci = cc_inv;
    for(int i = 0; i < d->number_of_correspondences; i++)
    {
        double * pi =  d->point_3d_array + 3 * i;
        double * a =  d->alphas_array + 4 * i;

        for(int j = 0; j < 3; j++)
        {
            a[1 + j] =
                    ci[3 * j    ] * (pi[0] - d->cws[0][0]) +
                    ci[3 * j + 1] * (pi[1] - d->cws[0][1]) +
                    ci[3 * j + 2] * (pi[2] - d->cws[0][2]);
            a[0] = 1.0 - a[1] - a[2] - a[3];
        }
    }
    // matrix_print(d->alphas_array,4,4,"Alphas array");


}
void  fill_M(epnp_data *d, double * M,
             const int row, const double * as, const double u, const double v)
{
    double * M1 = &M[0] + row * 12;
    double * M2 = M1 + 12;

    for(int i = 0; i < 4; i++) {
        M1[3 * i    ] = as[i] * d->focal_lenght_x;
        M1[3 * i + 1] = 0.0;
        M1[3 * i + 2] = as[i] * (d->center_x - u);

        M2[3 * i    ] = 0.0;
        M2[3 * i + 1] = as[i] * d->focal_lenght_y;
        M2[3 * i + 2] = as[i] * (d->center_y - v);
    }
    //// matrix_print(M,8,12,"M");
}



MAT_API void  getModelMat(double R[3][3], double t[3],float *out)
{
    float src[16];




    //assign t and r to ogl mat
    //rotation 1 row
    src[0]=(float)R[0][0];
    src[1]=(float)R[0][1];
    src[2]=(float)R[0][2];
    //translation 1 row
    src[3]=(float)t[0];


    src[4]=(float)R[1][0];
    src[5]=(float)R[1][1];
    src[6]=(float)R[1][2];
    src[7]=(float)t[1];


    src[8]= (float)R[2][0];
    src[9]= (float)R[2][1];
    src[10]=(float)R[2][2];
    src[11]=(float)t[2];
    src[12]=0.0f;
    src[13]=0.0f;
    src[14]=0.0f;
    src[15]=1.0f;
    //matrix_print_float( src,4,4,"++ROT++");


    float oglmat[16];
    matrix_zero_float(oglmat,16);
    oglmat[0*4]=1.0f;
    oglmat[1*4+1]=-1.0f;
    oglmat[2*4+2]=-1.0f;
    oglmat[3*4+3]=1.0f;

    matrix_mul_matrix_float(oglmat,src,out,4,4,4,4);

    matrix_transpose_float(out,4,4);

}




MAT_API double compute_pose(epnp_data *PnP,double R_est [3][3],double t_est[3])
{
    choose_control_points( PnP);
    compute_barycentric_coordinates(PnP);

    double M [2*4*12];

    for(int i = 0; i < 4; i++)
        fill_M( PnP ,M, 2 * i, PnP->alphas_array + 4 * i, PnP->point_2d_array[2 * i], PnP->point_2d_array[2 * i + 1]);



    // matrix_print(M,8,12,"MMMMMMMMMMMMMMMMMMMMMMMMM MAT");
    double Mt[12*8];


    matrix_transpose2(M,Mt,8,12);



    // matrix_print(Mt,12,8,getName(Mt));

    double resM[12*12];

    matrix_mul_matrix(Mt,M,resM,12,8,8,12);

    // matrix_print(resM,12,12,"A'*A");





    double U[12][12];
    double V[12][12];
    double singular_values[12];


    int error;//LIE IN THIS PLACE!!!
    error = Singular_Value_Decomposition(resM, 12, 12, (double*) U,
                                         singular_values, (double*) V);

    //for memory economy we assign to old unused variables
    matrix_transpose2(  (double*)U,resM,12,12);
     matrix_mul_scalar(resM,-1.0,12,12);

    // matrix_print(resM,12,12,"result after svd");
    // matrix_print(singular_values,1,12,getName(singular_values));
    //// matrix_print((double*)Ures,12,12,"result svd U");

    //cvMulTransposed(M, &MtM, 1);
    //cvSVD(&MtM, &D, &Ut, 0, CV_SVD_MODIFY_A | CV_SVD_U_T);
    // cvReleaseMat(&M);


    double L_6x10[6 * 10], Rho[6];
    matrix_zero(L_6x10,6*10);
    matrix_zero( Rho,6);
    //// part2
     // matrix_zero(resM,12*12);
     //matrix_load(resM,12*12);
    // matrix_print(resM,12,12,"LOADED M");

    compute_L_6x10(resM, L_6x10);
    // matrix_print(L_6x10,6,10,"l6x10");
    compute_rho(Rho,PnP);
    // matrix_print(Rho,1,6,"Rho");
    double Betas[4][4], rep_errors[4];
    double Rs[4][3][3], ts[4][3];

    find_betas_approx_1(L_6x10, Rho, Betas[1]);

    gapoint_2d_arrays_newton(L_6x10, Rho, Betas[1]);
    rep_errors[1] = compute_R_and_t(resM, Betas[1], Rs[1], ts[1],PnP);

    find_betas_approx_2(L_6x10, Rho, Betas[2]);
    gapoint_2d_arrays_newton(L_6x10, Rho, Betas[2]);
    rep_errors[2] = compute_R_and_t(resM, Betas[2], Rs[2], ts[2],PnP);

    find_betas_approx_3(L_6x10, Rho, Betas[3]);
    gapoint_2d_arrays_newton(L_6x10, Rho, Betas[3]);
    rep_errors[3] = compute_R_and_t(resM, Betas[3], Rs[3], ts[3],PnP);

    int N = 1;
    if (rep_errors[2] < rep_errors[1]) N = 2;
    if (rep_errors[3] < rep_errors[N]) N = 3;

    copy_R_and_t(Rs[N], ts[N], R_est, t_est);


    matrix_print((double*)R_est,3,3,"ROTATION");
    matrix_print((double*)t_est,1,3,"translation");


    return rep_errors[N];
}



