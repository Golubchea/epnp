#pragma once
#define MAT_API extern

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

MAT_API bool isRotationMatrix(float *R);
///
/// \brief rotationMatrixToEulerAngles
/// \param R[3*3] - matrix 3x3
/// \param angles[3] -Euler angles in degrees
///
MAT_API void rotationMatrixToEulerAngles_float(float * R,float *angles);// (float*)R[3][3] or R[3*3]


#ifdef __cplusplus
}
#endif // __cplusplus
