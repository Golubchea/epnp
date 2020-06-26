#include <pnp/epnp/epnp_data.h>
#include <pnp/epnp/epnp_rho_approx.h>
#include <pnp/epnp/math/singular_value_decomposition.h>
#include "math.h"
#include <pnp/epnp/math/matrix.h>


MAT_API void  copy_R_and_t(const double R_src[3][3], const double t_src[3],
double R_dst[3][3], double t_dst[3])
{
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++)
            R_dst[i][j] = R_src[i][j];
        t_dst[i] = t_src[i];
    }
}

MAT_API double  dist2(const double * p1, const double * p2)
{
    return
            (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]);
}

MAT_API double  dot(const double * v1, const double * v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


MAT_API void  estimate_R_and_t(double R[3][3], double t[3],epnp_data *d)
{
    double pc0[3], pw0[3];

    pc0[0] = pc0[1] = pc0[2] = 0.0;
    pw0[0] = pw0[1] = pw0[2] = 0.0;

    for(int i = 0; i < d->number_of_correspondences; i++) {
        const double * pc = d->pcs_array + 3 * i;
        const double * pw = d->point_3d_array + 3 * i;

        for(int j = 0; j < 3; j++) {
            pc0[j] += pc[j];
            pw0[j] += pw[j];
        }
    }
    for(int j = 0; j < 3; j++) {
        pc0[j] /= (double)d->number_of_correspondences;
        pw0[j] /= (double)d->number_of_correspondences;
    }

    double abt[3 * 3], abt_d[3], abt_u[3 * 3], abt_v[3 * 3];


    int var=0;
    for (var = 0; var < 9; ++var) {
        abt[var]=0;
    }
    for(int i = 0; i < d->number_of_correspondences; i++) {
        double * pc = d->pcs_array + 3 * i;
        double * pw = d->point_3d_array + 3 * i;

        for(int j = 0; j < 3; j++) {
            abt[3 * j    ] += (pc[j] - pc0[j]) * (pw[0] - pw0[0]);
            abt[3 * j + 1] += (pc[j] - pc0[j]) * (pw[1] - pw0[1]);
            abt[3 * j + 2] += (pc[j] - pc0[j]) * (pw[2] - pw0[2]);
        }
    }
    Singular_Value_Decomposition((double*)abt,3,3,(double*)abt_u,abt_d,(double*)abt_v);


    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            R[i][j] = dot(abt_u + 3 * i, abt_v + 3 * j);

    const double det =
            R[0][0] * R[1][1] * R[2][2] + R[0][1] * R[1][2] * R[2][0] + R[0][2] * R[1][0] * R[2][1] -
            R[0][2] * R[1][1] * R[2][0] - R[0][1] * R[1][0] * R[2][2] - R[0][0] * R[1][2] * R[2][1];

    if (det < 0) {
        R[2][0] = -R[2][0];
        R[2][1] = -R[2][1];
        R[2][2] = -R[2][2];
    }

    t[0] = pc0[0] - dot(R[0], pw0);
    t[1] = pc0[1] - dot(R[1], pw0);
    t[2] = pc0[2] - dot(R[2], pw0);
}

MAT_API void compute_ccs(const double * betas, const double * ut,epnp_data *d)
{
    for(int i = 0; i < 4; i++)
        d->ccs[i][0] = d->ccs[i][1] = d->ccs[i][2] = 0.0f;

    for(int i = 0; i < 4; i++) {
        const double * v = ut + 12 * (11 - i);
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 3; k++)
                d->ccs[j][k] += betas[i] * v[3 * j + k];
    }
}

MAT_API void compute_pcs(epnp_data *d)
{
    for(int i = 0; i < d->number_of_correspondences; i++) {
        double * a = d->alphas_array + 4 * i;
        double * pc = d->pcs_array + 3 * i;

        for(int j = 0; j < 3; j++)
            pc[j] = a[0] * d->ccs[0][j] + a[1] * d->ccs[1][j] + a[2] * d->ccs[2][j] + a[3] * d->ccs[3][j];
    }
}



MAT_API void  solve_for_sign(epnp_data *d)
{
    if (d->pcs_array[2] < 0.0) {
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < 3; j++)
                d->ccs[i][j] = -d->ccs[i][j];

        for(int i = 0; i < d->number_of_correspondences; i++) {
            d->pcs_array[3 * i    ] = -d->pcs_array[3 * i];
            d->pcs_array[3 * i + 1] = -d->pcs_array[3 * i + 1];
            d->pcs_array[3 * i + 2] = -d->pcs_array[3 * i + 2];
        }
    }
}


MAT_API double   reprojection_error(epnp_data *d, double R[3][3],  double t[3])
{
    double sum2 = 0.0;

    for(int i = 0; i < d->number_of_correspondences; i++) {
        double * pw = d->point_3d_array + 3 * i;
        double Xc = dot(R[0], pw) + t[0];
        double Yc = dot(R[1], pw) + t[1];
        double inv_Zc = 1.0 / (dot(R[2], pw) + t[2]);
        double ue = d->center_x + d->focal_lenght_x * Xc * inv_Zc;
        double ve = d->center_y + d->focal_lenght_y * Yc * inv_Zc;
        double u = d->point_2d_array[2 * i], v = d->point_2d_array[2 * i + 1];

        sum2 += sqrt( (u - ue) * (u - ue) + (v - ve) * (v - ve) );
    }

    return sum2 / d->number_of_correspondences;
}


MAT_API double  compute_R_and_t(const double * ut, const double * betas,
                                double R[3][3], double t[3],epnp_data *d)
{
    compute_ccs(betas, ut,d);
    compute_pcs(d);

    solve_for_sign(d);

    estimate_R_and_t(R, t,d);

    return  reprojection_error(d,R, t);
}






// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_1 = [B11 B12     B13         B14]

MAT_API void  find_betas_approx_1(const double * L_6x10, const double * Rho,
                                  double * betas)
{
    double l_6x4[6 * 4], b4[4];


    for(int i = 0; i < 6; i++)
    {
        l_6x4[4*i+0]=L_6x10[10*i+0];
        l_6x4[4*i+1]=L_6x10[10*i+1];
        l_6x4[4*i+2]=L_6x10[10*i+3];
        l_6x4[4*i+3]=L_6x10[10*i+6];

    }


    double U[6][4];
    double V[6][4];
    double D[4];



    Singular_Value_Decomposition((double*) l_6x4, 6, 4, (double*) U,  D, (double*) V);

 //   matrix_print((double*)U,6,4,"find betas 1 svd U\n");
 //   matrix_print((double*)V,6,4,"find betas 1 svd V\n");
 //   matrix_print((double*)D,1,4,"find betas 1 svd D\n");
    double tolerance=0.000001;
    Singular_Value_Decomposition_Solve((double*)U,D,(double*)V,tolerance,6,4,Rho,b4)  ;

  //  matrix_print( b4,1,4,"find betas 1 svd solve b4\n");



    if (b4[0] < 0) {
        betas[0] = sqrt(-b4[0]);
        betas[1] = -b4[1] / betas[0];
        betas[2] = -b4[2] / betas[0];
        betas[3] = -b4[3] / betas[0];
    } else {
        betas[0] = sqrt(b4[0]);
        betas[1] = b4[1] / betas[0];
        betas[2] = b4[2] / betas[0];
        betas[3] = b4[3] / betas[0];
    }

  //  matrix_print((double*)betas,1,4,"find betas 1 svd solve b4\n");
}

// betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
// betas_approx_2 = [B11 B12 B22                            ]

MAT_API void  find_betas_approx_2(const double * L_6x10, const double * Rho,
                                  double * betas)
{
    double l_6x3[6 * 3], b3[3];


    for(int i = 0; i < 6; i++)
    {
        l_6x3[3*i+0]=L_6x10[10*i+0];
        l_6x3[3*i+1]=L_6x10[10*i+1];
        l_6x3[3*i+2]=L_6x10[10*i+2];

    }

    double U[6][3];
    double V[6][3];
    double D[3];



    Singular_Value_Decomposition((double*) l_6x3, 6, 3, (double*) U,  D, (double*) V);
    double tolerance=0.000001;
    Singular_Value_Decomposition_Solve((double*)U,D,(double*)V,tolerance,6,3,Rho,b3)  ;


    if (b3[0] < 0) {
        betas[0] = sqrt(-b3[0]);
        betas[1] = (b3[2] < 0) ? sqrt(-b3[2]) : 0.0;
    } else {
        betas[0] = sqrt(b3[0]);
        betas[1] = (b3[2] > 0) ? sqrt(b3[2]) : 0.0;
    }

    if (b3[1] < 0) betas[0] = -betas[0];

    betas[2] = 0.0;
    betas[3] = 0.0;
}


MAT_API void  find_betas_approx_3(const double * L_6x10, const double * Rho,
                                  double * betas)
{
    double l_6x5[6 * 5], b5[5];


    for(int i = 0; i < 6; i++)
    {
        l_6x5[5*i+0]=L_6x10[10*i+0];
        l_6x5[5*i+1]=L_6x10[10*i+1];
        l_6x5[5*i+2]=L_6x10[10*i+2];
        l_6x5[5*i+3]=L_6x10[10*i+3];
        l_6x5[5*i+4]=L_6x10[10*i+4];

    }



    double U[6][5];
    double V[6][5];
    double D[5];



    Singular_Value_Decomposition((double*) l_6x5, 6, 5, (double*) U,  D, (double*) V);
    double tolerance=0.000001;
    Singular_Value_Decomposition_Solve((double*)U,D,(double*)V,tolerance,6,5,Rho,b5)  ;

    if (b5[0] < 0) {
        betas[0] = sqrt(-b5[0]);
        betas[1] = (b5[2] < 0) ? sqrt(-b5[2]) : 0.0;
    } else {
        betas[0] = sqrt(b5[0]);
        betas[1] = (b5[2] > 0) ? sqrt(b5[2]) : 0.0;
    }
    if (b5[1] < 0) betas[0] = -betas[0];
    betas[2] = b5[3] / betas[0];
    betas[3] = 0.0;
}

MAT_API void  compute_L_6x10(const double * ut, double * l_6x10)
{
    const double * v[4];

    v[0] = ut + 12 * 11;
    v[1] = ut + 12 * 10;
    v[2] = ut + 12 *  9;
    v[3] = ut + 12 *  8;

    double dv[4][6][3];

    for(int i = 0; i < 4; i++) {
        int a = 0, b = 1;
        for(int j = 0; j < 6; j++) {
            dv[i][j][0] = v[i][3 * a    ] - v[i][3 * b];
            dv[i][j][1] = v[i][3 * a + 1] - v[i][3 * b + 1];
            dv[i][j][2] = v[i][3 * a + 2] - v[i][3 * b + 2];

            b++;
            if (b > 3) {
                a++;
                b = a + 1;
            }
        }
    }

    for(int i = 0; i < 6; i++) {
        double * row = l_6x10 + 10 * i;

        row[0] =        dot(dv[0][i], dv[0][i]);
        row[1] = 2.0 * dot(dv[0][i], dv[1][i]);
        row[2] =        dot(dv[1][i], dv[1][i]);
        row[3] = 2.0 * dot(dv[0][i], dv[2][i]);
        row[4] = 2.0 * dot(dv[1][i], dv[2][i]);
        row[5] =        dot(dv[2][i], dv[2][i]);
        row[6] = 2.0 * dot(dv[0][i], dv[3][i]);
        row[7] = 2.0 * dot(dv[1][i], dv[3][i]);
        row[8] = 2.0 * dot(dv[2][i], dv[3][i]);
        row[9] =        dot(dv[3][i], dv[3][i]);
    }
}

MAT_API void  compute_rho(double * rho ,epnp_data *d)
{
    rho[0] = dist2(d->cws[0], d->cws[1]);
    rho[1] = dist2(d->cws[0], d->cws[2]);
    rho[2] = dist2(d->cws[0], d->cws[3]);
    rho[3] = dist2(d->cws[1], d->cws[2]);
    rho[4] = dist2(d->cws[1], d->cws[3]);
    rho[5] = dist2(d->cws[2], d->cws[3]);
}

MAT_API void  compute_A_and_b_gapoint_2d_arrays_newton(const double * l_6x10, const double * rho,
                                                       double betas[4], double * A, double * b)
{
    for(int i = 0; i < 6; i++)
    {
        const double * rowL = l_6x10 + i * 10;
        double * rowA =A + i * 4;




        rowA[0] = 2.0 * rowL[0] * betas[0] +     rowL[1] * betas[1] +     rowL[3] * betas[2] +     rowL[6] * betas[3];
        rowA[1] =     rowL[1] * betas[0] + 2.0 * rowL[2] * betas[1] +     rowL[4] * betas[2] +     rowL[7] * betas[3];
        rowA[2] =     rowL[3] * betas[0] +     rowL[4] * betas[1] + 2.0 * rowL[5] * betas[2] +     rowL[8] * betas[3];
        rowA[3] =     rowL[6] * betas[0] +     rowL[7] * betas[1] +     rowL[8] * betas[2] + 2.0 * rowL[9] * betas[3];



        b[i]=rho[i] -
                (
                    rowL[0] * betas[0] * betas[0] +
                rowL[1] * betas[0] * betas[1] +
                rowL[2] * betas[1] * betas[1] +
                rowL[3] * betas[0] * betas[2] +
                rowL[4] * betas[1] * betas[2] +
                rowL[5] * betas[2] * betas[2] +
                rowL[6] * betas[0] * betas[3] +
                rowL[7] * betas[1] * betas[3] +
                rowL[8] * betas[2] * betas[3] +
                rowL[9] * betas[3] * betas[3]
                );

    }
}


MAT_API void  qr_solve(double * A, double * b, double * X,int m,int n)
{

//    matrix_print(A,6,4,"aaa");
//    matrix_print(b,1,6,"b");
//    matrix_print(b,1,6,"b");


    static int max_nr = 0;
    static double * A1 , * A2 ;

    const int nr = m;
    const int nc =n;

    if (max_nr != 0 && max_nr < nr)
    {
        free(A1) ;
        free(A2);
    }
    if (max_nr < nr)
    {
        max_nr = nr;
        A1 = (double*)malloc(nr*sizeof (double)  );
        matrix_zero(A1,nr);
        A2 = (double*)malloc(nr*sizeof (double)  );
        matrix_zero(A2,nr);
    }

    double * pA = &A[0], * ppAkk = pA;
    for(int k = 0; k < nc; k++) {
        double * ppAik = ppAkk, eta = fabs(*ppAik);
        for(int i = k + 1; i < nr; i++) {
            double elt = fabs(*ppAik);
            if (eta < elt) eta = elt;
            ppAik += nc;
        }

        if (eta == 0) {
            A1[k] = A2[k] = 0.0;
            //printf("God damnit, A is singular, this shouldn't happen. \n");
            return;
        } else {
            double * ppAik = ppAkk, sum = 0.0, inv_eta = 1. / eta;
            for(int i = k; i < nr; i++) {
                *ppAik *= inv_eta;
                sum += *ppAik * *ppAik;
                ppAik += nc;
            }
            double sigma = sqrt(sum);
            if (*ppAkk < 0)
                sigma = -sigma;
            *ppAkk += sigma;
            A1[k] = sigma * *ppAkk;
            A2[k] = -eta * sigma;
            for(int j = k + 1; j < nc; j++) {
                double * ppAik = ppAkk, sum = 0;
                for(int i = k; i < nr; i++) {
                    sum += *ppAik * ppAik[j - k];
                    ppAik += nc;
                }
                double tau = sum / A1[k];
                ppAik = ppAkk;
                for(int i = k; i < nr; i++) {
                    ppAik[j - k] -= tau * *ppAik;
                    ppAik += nc;
                }
            }
        }
        ppAkk += nc + 1;
    }

    // b <- Qt b
    double * ppAjj = pA, * pb = &b[0];
    for(int j = 0; j < nc; j++) {
        double * ppAij = ppAjj, tau = 0;
        for(int i = j; i < nr; i++)	{
            tau += *ppAij * pb[i];
            ppAij += nc;
        }
        tau /= A1[j];
        ppAij = ppAjj;
        for(int i = j; i < nr; i++) {
            pb[i] -= tau * *ppAij;
            ppAij += nc;
        }
        ppAjj += nc + 1;
    }

    // X = R-1 b
    double * pX = X;
    pX[nc - 1] = pb[nc - 1] / A2[nc - 1];
    for(int i = nc - 2; i >= 0; i--) {
        double * ppAij = pA + i * nc + (i + 1), sum = 0;

        for(int j = i + 1; j < nc; j++) {
            sum += *ppAij * pX[j];
            ppAij++;
        }
        pX[i] = (pb[i] - sum) / A2[i];
    }


}

MAT_API void  gapoint_2d_arrays_newton(const double * L_6x10, const double * Rho,
                                       double betas[4])
{
    const int iterations_number = 5;

    double a[6*4], b[6], x[4];
    matrix_zero(a,6*4);
    matrix_zero(b,6);
    matrix_zero(x,4);


    for(int k = 0; k < iterations_number; k++)
    {
        compute_A_and_b_gapoint_2d_arrays_newton(L_6x10, Rho,
                                                 betas, a, b);


        qr_solve(a, b, x,6,4);

        for(int i = 0; i < 4; i++)
            betas[i] += x[i];
    }
}


