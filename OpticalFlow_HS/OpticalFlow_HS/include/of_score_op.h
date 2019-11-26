/* ------------------------------ */
/* --- OpticalFlow_operator.h --- */
/* ------------------------------ */

/*
 * Copyright (c) 2017, Lionel Lacassagne, all rights reserved
 * LIP6, UPMC, CNRS
 */

#ifndef _OPTICAL_FLOW_OPERATOR_H_
#define _OPTICAL_FLOW_OPERATOR_H_

#ifdef __cplusplus
extern "C" {
#endif


float32 calcMSE_F32(float32 **X, int height, int width, float32 **Y);
float32 calcRMSE_F32(float32 **X, int height, int width, float32 **Y);

float32 calcMSE_I8 (uint8    **X, int height, int width, uint8    **Y);
double calcAngle(double x0, double y0, double x1, double y1);
float32 calcAngle_F32(float32 **X0, float32 **Y, int height, int width, float32 **X1, float32 **Y1);

// convert deegre into radian
double of_rad2deg(double rad);

// compute AE between two vectors (u0,v0) and (u1,v1)
double of_angularError_F32(float32 u0, float32 v0, float32 u1, float32 v1);

// compute AE between two vectors (u0,v0) and (u1,v1)
double of_endpointError_F32(float32 u0, float32 v0, float32 u1, float32 v1);

// compute AE between two matrices of vectors (U0,V0) and (U1,V1)
double of_calcAngularError_F32(float32 **U0, float32 **V0, float32 **U1, float32 **V1, int height, int width);

// compute EE between two matrices of vectors (U0,V0) and (U1,V1)
double of_calcEndpointError_F32(float32 **U0, float32 **V0, float32 **U1, float32 **V1, int height, int width);

// compute IE between input image and ground-true / reconstructed image
// compute IE = RMS Error between input image and ground-true / reconstructed image
double of_calcSquareError_F32        (float32 **I0, float32 **I1, int height, int width);
double of_calcMeanSquareError_F32    (float32 **I0, float32 **I1, int height, int width);
double of_calcMeanSquareError_U16    (uint16  **I0, uint16  **I1, int height, int width);
double of_calcRootMeanSquareError_F32(float32 **I0, float32 **I1, int height, int width);
double of_calcInterpolationError_F32 (float32 **I0, float32 **I1, int height, int width);

double of_calcMeanSquareError_U16_F32(uint16 **I0, int qi, float32 **I1, int height, int width);
double of_calcMeanSquareError_S16_F32(sint16 **I0, int qi, float32 **I1, int height, int width);
double of_calcMeanSquareError_S32_F32(sint32 **I0, int qi, float32 **I1, int height, int width);

double og_calcRelativeError_F32(float32 **Xref, float32 **X, int height, int width);
#ifdef __cplusplus
}
#endif

#endif // _OPTICAL_FLOW_OPERATOR_H_