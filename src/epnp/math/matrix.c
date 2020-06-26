
#include <pnp/epnp/math/matrix.h>
#include <stdbool.h>
MAT_API inline void matrix_mul_scalar(double * src ,double scalar,unsigned int M,unsigned int N)
{
    unsigned i=0;


    for (i = 0; i <M*N; ++i)
    {
        src[i]=src[i]*scalar;
    }


}


MAT_API inline void matrix_transpose2(double  *array,double  *out , int m, int n)
{

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int old_idx = i * n + j;
            int new_idx = j * m + i;
            out[new_idx] = array[old_idx];
        }
    }

}

MAT_API inline void matrix_transpose(double  *array, int m, int n)
{
    double new_array[m * n];
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int old_idx = i * n + j;
            int new_idx = j * m + i;
            new_array[new_idx] = array[old_idx];
        }
    }
    for (int i = 0; i < m * n; i++)
    {
        array[i] = new_array[i];
    }
}



MAT_API inline void matrix_transpose_float(float  *array, int m, int n)
{
    float new_array[m * n];
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            int old_idx = i * n + j;
            int new_idx = j * m + i;
            new_array[new_idx] = array[old_idx];
        }
    }
    for (int i = 0; i < m * n; i++)
    {
        array[i] = new_array[i];
    }
}


//   A[4*0+0] ; A[4*0+1] ; A[4*0+2] ; A[4*0+3];
//   A[4*1+0] ; A[4*1+1] ; A[4*1+2] ; A[4*1+3];
//   A[4*2+0] ; A[4*2+1] ; A[4*2+2] ; A[4*2+3];

//    B[3*0+0]=0.5;     B[3*0+1]=0;     B[3*0+2]=-0.5;
//    B[3*1+0]=0.5;     B[3*1+1]=0;     B[3*1+2]=0.5;
//    B[3*2+0]=-0.5;    B[3*2+1]=0;     B[3*2+2]=-0.5;
//    B[3*3+0]=-0.5;    B[3*3+1]=0;     B[3*3+2]=0.5;


//    Arows,  Acols,  Brows,  Bcols
//mat 3        4    x  4        3      = 3x3
//Arows Bcols


//    res[ Bcols*0+0]=A[Acols*0+0]*B[Bcols*0+0]+A[Acols*0+1]*B[Bcols*1+0]+A[Acols*0+2]*B[Bcols*2+0]+A[Acols*0+3]*B[Bcols*3+0];
//    res[ Bcols*0+1]=A[Acols*0+0]*B[Bcols*0+1]+A[Acols*0+1]*B[Bcols*1+1]+A[Acols*0+2]*B[Bcols*2+1]+A[Acols*0+3]*B[Bcols*3+1];
//    res[ Bcols*0+2]=A[Acols*0+0]*B[Bcols*0+2]+A[Acols*0+1]*B[Bcols*1+2]+A[Acols*0+2]*B[Bcols*2+2]+A[Acols*0+3]*B[Bcols*3+2];

//    res[ Bcols*1+0]=A[Acols*1+0]*B[Bcols*0+0]+A[Acols*1+1]*B[Bcols*1+0]+A[Acols*1+2]*B[Bcols*2+0]+A[Acols*1+3]*B[Bcols*3+0];
//    res[ Bcols*1+1]=A[Acols*1+0]*B[Bcols*0+1]+A[Acols*1+1]*B[Bcols*1+1]+A[Acols*1+2]*B[Bcols*2+1]+A[Acols*1+3]*B[Bcols*3+1];
//    res[ Bcols*1+2]=A[Acols*1+0]*B[Bcols*0+2]+A[Acols*1+1]*B[Bcols*1+2]+A[Acols*1+2]*B[Bcols*2+2]+A[Acols*1+3]*B[Bcols*3+2];

//    res[ Bcols*2+0]=A[Acols*2+0]*B[Bcols*0+0]+A[Acols*2+1]*B[Bcols*1+0]+A[Acols*2+2]*B[Bcols*2+0]+A[Acols*2+3]*B[Bcols*3+0];
//    res[ Bcols*2+1]=A[Acols*2+0]*B[Bcols*0+1]+A[Acols*2+1]*B[Bcols*1+1]+A[Acols*2+2]*B[Bcols*2+1]+A[Acols*2+3]*B[Bcols*3+1];
//    res[ Bcols*2+2]=A[Acols*2+0]*B[Bcols*0+2]+A[Acols*2+1]*B[Bcols*1+2]+A[Acols*2+2]*B[Bcols*2+2]+A[Acols*2+3]*B[Bcols*3+2];



MAT_API inline void matrix_mul_matrix(double * A ,double * B,double *res,unsigned int Arows,unsigned int Acols,unsigned int Brows,unsigned int Bcols)
{
    assert(Acols==Brows);

    unsigned i=0;
    unsigned j=0;
    unsigned k=0;
    for ( i = 0; i <  Arows; i++)
    {
        for ( j = 0; j < Bcols; j++)
        {
            double sum = 0.0;
            for (k = 0; k < Brows; k++)
            {
                sum = sum + A[i * Acols + k] * B[k * Bcols + j];
            }
            res[i * Bcols + j] = sum;
        }
    }
}







MAT_API inline void matrix_add_matrix(double * A ,double * B,double *res,unsigned int Arows,unsigned int Acols)
{


    unsigned i=0;
    unsigned j=0;
    unsigned k=0;
    for ( i = 0; i <  Arows; i++)
    {
        for ( j = 0; j < Acols; j++)
        {

            res[i *Arows + j] = A[i *Arows + j]  +B[i *Arows + j];
        }
    }
}


MAT_API inline void matrix_substract_matrix(double * A ,double * B,double *res,unsigned int Arows,unsigned int Acols)
{


    unsigned i=0;
    unsigned j=0;
    unsigned k=0;
    for ( i = 0; i <  Arows; i++)
    {
        for ( j = 0; j < Acols; j++)
        {

            res[i *Arows + j] = A[i *Arows + j]  -  B[i *Arows + j];
        }
    }
}


MAT_API inline void matrix_mul_matrix_float(float * A ,float * B,float *res,unsigned int Arows,unsigned int Acols,unsigned int Brows,unsigned int Bcols)
{
    assert(Acols==Brows);

    unsigned i=0;
    unsigned j=0;
    unsigned k=0;
    for ( i = 0; i <  Arows; i++)
    {
        for ( j = 0; j < Bcols; j++)
        {
            float sum = 0.0;
            for (k = 0; k < Brows; k++)
            {
                sum = sum + A[i * Acols + k] * B[k * Bcols + j];
            }
            res[i * Bcols + j] = sum;
        }
    }
}

MAT_API void Atranspose_mul_A(double *src,double *res,unsigned  M,unsigned   N) //use malloc res[N+1]
{
    double *A_transposed=(double*)malloc((N )*(M)*sizeof (double));

    matrix_transpose2(src,A_transposed,M,N);
    matrix_mul_matrix(A_transposed,src,res,N,M,M,N);
    free(  A_transposed);
}

MAT_API void matrix_print(double* c,unsigned int M, unsigned int N, char* name)
{
    printf("%s matrix%ux%u :\n",name,M,N);
    for(unsigned int i = 0; i <M; i++)
    {
        for(unsigned int j = 0; j < N; j++)
        {
            printf("%lf ",c[i*N+j]);
        }
        printf("\n");
    }
}

MAT_API void matrix_print_float(float* c,unsigned int M, unsigned int N, char* name)
{
    printf("%s matrix%ux%u :\n",name,M,N);
    for(unsigned int i = 0; i <M; i++)
    {
        for(unsigned int j = 0; j < N; j++)
        {
            printf("%g ",c[i*N+j]);
        }
        printf("\n");
    }
}

MAT_API void matrix_load(double * mat, int M)
{
    FILE *f;
    f = fopen("mat.txt", "r");

    for (int i = 0; i < M; ++i)
    {
        fscanf(f, "%lf;", &mat[i]);
    }

}
MAT_API void matrix_zero(double * mat, int M)
{

    for (int i = 0; i < M; ++i)
    {
        mat[i]=0.0;
    }

}

MAT_API void matrix_copy(double * src,double * dst, int rang)
{

    for (int i = 0; i < rang*rang; ++i)
    {
        dst[i]=src[i];
    }

}


MAT_API void matrix_copy_float(float * src,float * dst, int rang)
{

    for (int i = 0; i < rang*rang; ++i)
    {
        dst[i]=src[i];
    }

}




MAT_API void matrix_identity(double * mat, int M)
{
    matrix_zero(mat,M*M);

    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            mat[i*M+i]=1.0 ;
        }

    }
}

MAT_API void matrix_zero_float(float * mat, int M)
{
    for (int i = 0; i < M; ++i)
    {
        mat[i]=0.0;
    }
}

MAT_API void matrix_identity_float(float * mat, int M)
{
    matrix_zero_float(mat,M*M);

    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            mat[i*M+i]=1.0f;
        }

    }
}


MAT_API void matrix_inverse3x3(double* m,double * res)
{
    double det = m[0] * (m[1*3+ 1] * m[3*2+ 2] - m[3*2+ 1] * m[3*1+ 2]) -
            m[0+ 1] * (m[1*3+ 0] * m[2*3+ 2] - m[3*1+ 2] * m[2*3+ 0]) +
            m[0+ 2] * (m[1*3+ 0] * m[2*3+ 1] - m[1*3+ 1] * m[2*3+ 0]);

    double invdet = 1.0 / det;

    printf("det %f \n",det);
    res[0] = (m[1*3+ 1] * m[2*3+ 2] - m[2*3+ 1] * m[1*3+ 2]) * invdet;
    res[1] = (m[0*3+ 2] * m[2*3+ 1] - m[0*3+ 1] * m[2*3+ 2]) * invdet;
    res[2] = (m[0*3+ 1] * m[1*3+ 2] - m[0*3+ 2] * m[1*3+ 1]) * invdet;
    res[3] = (m[1*3+ 2] * m[2*3+ 0] - m[1*3+ 0] * m[2*3+ 2]) * invdet;
    res[4] = (m[0*3+ 0] * m[2*3+ 2] - m[0*3+ 2] * m[2*3+ 0]) * invdet;
    res[5] = (m[1*3+ 0] * m[0*3+ 2] - m[0*3+ 0] * m[1*3+ 2]) * invdet;
    res[6] = (m[1*3+ 0] * m[2*3+ 1] - m[2*3+ 0] * m[1*3+ 1]) * invdet;
    res[7] = (m[2*3+ 0] * m[0*3+ 1] - m[0*3+ 0] * m[2*3+ 1]) * invdet;
    res[8] = (m[0*3+ 0] * m[1*3+ 1] - m[1*3+ 0] * m[0*3+ 1]) * invdet;
}






MAT_API bool matrix_inverse4x4(const double * m , double *invOut )
{
    double inv[16], det;


    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;
    int i;
    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}


MAT_API bool matrix_inverse4x4_float( float *m)
{
    float inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0f / det;


    matrix_zero_float(m,4*4);
    for (i = 0; i < 16; i++)
       m[i] = inv[i] * det;

    return true;
}



MAT_API inline void matrix_translate4x4_float(float * srcmat4x4 ,float * translate_vector3)
{
    float translate_mat[16];

    matrix_identity_float(translate_mat,4);

    translate_mat[12]=translate_vector3[0];
    translate_mat[13]=translate_vector3[1];
    translate_mat[14]=translate_vector3[2];

    float src_temp[16];
    matrix_copy_float( srcmat4x4,src_temp,4);

    // matrix_print_float(translate_mat,4,4,"floatmat");
    matrix_mul_matrix_float(translate_mat,src_temp,srcmat4x4,4,4,4,4);

}


