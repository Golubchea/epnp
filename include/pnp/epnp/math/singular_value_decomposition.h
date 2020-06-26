////////////////////////////////////////////////////////////////////////////////
// File: singular_value_decomposition.c                                       //
// Contents:                                                                  //
//    Singular_Value_Decomposition                                            //
//    Singular_Value_Decomposition_Solve                                      //
//    Singular_Value_Decomposition_Inverse                                    //
////////////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

int Singular_Value_Decomposition(double* A, int nrows, int ncols, double* U,
                                 double* singular_values, double* V);



void Singular_Value_Decomposition_Solve(double* U, double* D, double* V,
                                        double tolerance, int nrows, int ncols, double *B, double* x) ;

void Singular_Value_Decomposition_Inverse(double* U, double* D, double* V,
                                          double tolerance, int nrows, int ncols, double *Astar) ;

#ifdef __cplusplus
}
#endif // __cplusplus
                                                                          //



