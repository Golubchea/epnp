
#pragma once
#define MAT_API extern
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus
//all matMxN[i][j] in form matMxN[i*M+j]

// check this function
//+
MAT_API inline void matrix_add_matrix(double * A ,double * B,double *res,unsigned int Arows,unsigned int Acols);
//-
MAT_API inline void matrix_substract_matrix(double * A ,double * B,double *res,unsigned int Arows,unsigned int Acols);

// can parallelize this functions
MAT_API inline void matrix_mul_matrix(double * A ,double * B,double *res,unsigned int Arows,unsigned int Acols,unsigned int Brows,unsigned int Bcols);
MAT_API inline void matrix_mul_matrix_float(float * A ,float * B,float *res,unsigned int Arows,unsigned int Acols,unsigned int Brows,unsigned int Bcols);

MAT_API inline void matrix_mul_scalar(double * src ,double scalar,unsigned int M,unsigned int N);

MAT_API inline void matrix_transpose(double  *array, int m, int n);
MAT_API inline void matrix_transpose_float(float  *array, int m, int n);
MAT_API inline void matrix_transpose2(double  *array,double  *out , int m, int n);
//use malloc res[N+1]
MAT_API void Atranspose_mul_A(double *src,double *res,unsigned  M,unsigned   N);

//debug print
MAT_API void matrix_print(double* c,unsigned int M, unsigned int N, char* name);
MAT_API void matrix_print_float(float* c,unsigned int M, unsigned int N, char* name);
MAT_API void matrix_load(double * mat, int M);
MAT_API void matrix_zero(double * mat, int M);
MAT_API void matrix_zero_float(float * mat, int M);


MAT_API void matrix_copy(double * src,double * dst, int rang);
MAT_API void matrix_copy_float(float * src,float * dst, int rang);

MAT_API void matrix_identity(double * mat, int M);
MAT_API void matrix_identity_float(float * mat, int M);

MAT_API void matrix_inverse3x3(double* m,double * res);
MAT_API bool matrix_inverse4x4(const double *m , double *invOut);


MAT_API bool matrix_inverse4x4_float( float *m);


//translate object
MAT_API inline void matrix_translate4x4_float(float * srcmat4x4 ,float * translate_vector3);


#ifdef __cplusplus
}
#endif // __cplusplus





