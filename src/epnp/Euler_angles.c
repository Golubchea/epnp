//R is 3x3 array
#include <pnp/epnp/math/matrix.h>
#include <stdbool.h>
#include <math.h>
bool isRotationMatrix(float *R)
{
    float Rt[3*3];
    matrix_transpose2(R,Rt,3,3);

    float shouldBeIdentity[3*3];
    matrix_mul_matrix(Rt,R,shouldBeIdentity,3,3,3,3);

    float I[3*3]={1,0,0,   0,1,0,   0,0,1};


    float res[3*3];
    // R^t*R-E=0
     matrix_substract_matrix(shouldBeIdentity,I,res,3,3);
     //TODO : make norm
     //norm(I, shouldBeIdentity) < 1e-6;

     if(res[0]<1e-6)
     {
         return  1;
     }
     else
     {
         return 0;
     }

}

MAT_API void rotationMatrixToEulerAngles_float(float * R,float *angles)// (float*)R[4][4] or R[4*4]
{
 
   // assert(isRotationMatrix(R));
     float Yaw,Pitch,Roll;
    if (R[0*4+0] == 1.0f)
            {
                Yaw = atan2f(R[0*4+2], R[2*4+3] );
                Pitch = 0;
                Roll = 0;

            }else if (R[0*4+0] == -1.0f)
            {
                Yaw = atan2f(R[0*4 +2], R[2*4+3]);
                Pitch = 0;
                Roll = 0;
            }else
            {

                Yaw = atan2f(-R[2*4+0],R[0*4+0]);
                Pitch = asinf(R [1*4+0]);
                Roll = atan2f(-R[1*4+2],R[1*4+1]);
            }


    angles[0]=Yaw/(float)M_PI *180.0f;
    angles[1]=Pitch/(float)M_PI *180.0f;
    angles[2]=Roll/(float)M_PI *180.0f;
}
