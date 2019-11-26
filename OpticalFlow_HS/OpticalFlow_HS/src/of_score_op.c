/* ------------------------------------ */
/* --- OpticalFlow_score_operator.c --- */
/* ------------------------------------ */

/*
 * Copyright (c) 2017, Lionel Lacassagne, all rights reserved
 * LIP6, UPMC, CNRS
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// NRC lib
#include "nrc.h"
#include "nrdef.h"
#include "nrtype.h"
#include "nralloc.h"
#include "nrarith.h"
#include "nrio.h"

// local
#include "of_score_op.h"

// --------------------------
double of_rad2deg(double rad)
// --------------------------
{
    // convert deegre into radian
    
    return 180.0 * rad / 3.141592653589793;
}
/* -------------------------------------------------------------- */
float32 calcMSE_F32(float32 **X, int height, int width, float32 **Y)
/* -------------------------------------------------------------- */
{
    int i, j;
    
    float32 x, y;
    float32 d, e = 0.0f;
    float32 mse, sz;

    sz = (float32) height * width;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) {
            x = X[i][j];
            y = Y[i][j];
            
            d = (x - y);
            e = e + d * d;
        }
    }
    //mse = (float32) sqrt(e) / sz;
    mse = (float32) e / sz;
    return mse;
}
/* -------------------------------------------------------------- */
float32 calcRMSE_F32(float32 **X, int height, int width, float32 **Y)
/* -------------------------------------------------------------- */
{
    float32 mse = calcMSE_F32(X, height, width, Y);
    float32 rmse = sqrtf(mse);
    return rmse;
}
/* ------------------------------------------------------- */
float32 calcMSE_I8(uint8 **X, int height, int width, uint8 **Y)
/* ------------------------------------------------------- */
{
    int i, j;
    
    float32 x, y;
    float32 d, e = 0.0f;
    float32 mse, sz;
    
    sz = (float32) height * width;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) {
            x = (float32) X[i][j];
            y = (float32) Y[i][j];
            d = (x - y);
            e = e + d * d;
        }
    }
    mse = (float32) sqrt(e) / sz;
    return mse;
}
/* ------------------------------------------------------ */
double calcAngle(double x0, double y0, double x1, double y1)
/* ------------------------------------------------------ */
{
    double n2;
    double dot;
    double dot2;
    double c, c2;
    double ac;
    
    ac = 0.0;
    dot = x0 * x1 + y0 * y1;
    dot2 = dot * dot;
    n2 = (x0*x0 + y0*y0) * (x1*x1 + y1*y1);
    
    
    if(n2 != 0.0) {
        
        c2 = dot2 / n2;
        c  = sqrt(c2);
        if(dot < 0.0) c = -c;
        ac = acos(c);
    }
    ac = (180 * ac) / PI;
    return ac;
}
/* ---------------------------------------------------------------------------------------------- */
float32 calcAngle_F32(float32 **X0, float32 **Y0, int height, int width, float32 **X1, float32 **Y1)
/* ---------------------------------------------------------------------------------------------- */
{
    int i, j;
    
    float64 x0, y0;
    float64 x1, y1;
    float64 a, sa, sad;
    float64 sz;
    
    sz = (float32) height * width;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) {
            x0 = X0[i][j];
            y0 = Y0[i][j];
            
            x1 = X1[i][j];
            y1 = Y1[i][j];
            
            a = calcAngle(x0, y0, x1, y1);
            
            sa += a;
        }
    }
    sad = (float32) (sa / sz);
    
    return sad;
}
// -----------------------------------------------------------------------
double of_angularError_F32(float32 u0, float32 v0, float32 u1, float32 v1)
// -----------------------------------------------------------------------
{
    // compute AE between two vectors (u0,v0) and (u1,v1)
    
    double num = 1.0 + u0 * u1 + v0 * v1;
    double den = sqrt((1.0 + u0 * u0 + v0 * v0) * (1.0 + u1 * u1 + v1 * v1));
    double rad = acos(num / den);
    double deg = of_rad2deg(rad);
    return deg;
}
// ------------------------------------------------------------------------
double of_endpointError_F32(float32 u0, float32 v0, float32 u1, float32 v1)
// ------------------------------------------------------------------------
{
    // compute AE between two vectors (u0,v0) and (u1,v1)
    
    double d = (u0 - u1) * (u0 - u1) + (v0 - v1) * (v0 - v1);
    double ee = sqrt(d);
    return ee;
}
// ----------------------------------------------------------------------------------------------------------
double of_calcAngularError_F32(float32 **U0, float32 **V0, float32 **U1, float32 **V1, int height, int width)
// ----------------------------------------------------------------------------------------------------------
{
    // compute AE between two matrices of vectors (U0,V0) and (U1,V1)
    
    double e, sx, s; // pixel, line, total error
    
    s = 0.0;
    for(int i = 0; i < height; i++) {
        sx = 0.0;
        for(int j = 0; j< width; j++) {
            e = of_angularError_F32(U0[i][j], V0[i][j], U1[i][j], V1[i][j]);
            sx += e;
        }
        s += sx;
    }
    s = s / (height * width);
    return s;
}
// -----------------------------------------------------------------------------------------------------------
double of_calcEndpointError_F32(float32 **U0, float32 **V0, float32 **U1, float32 **V1, int height, int width)
// -----------------------------------------------------------------------------------------------------------
{
    // compute EE between two matrices of vectors (U0,V0) and (U1,V1)
    
    double e, sx, s; // pixel, line, total error
    
    s = 0.0;
    for(int i = 0; i < height; i++) {
        sx = 0.0;
        for(int j = 0; j< width; j++) {
            e = of_endpointError_F32(U0[i][j], V0[i][j], U1[i][j], V1[i][j]);
            sx += e;
        }
        s += sx;
    }
    s = s / (height * width);
    return s;
}
// -----------------------------------------------------------------------------
double of_calcSquareError_F32(float32 **I0, float32 **I1, int height, int width)
// -----------------------------------------------------------------------------
{
    double e, sx, s; // pixel, line, total error
    
    s = 0.0;
    
    for(int i = 0; i < height; i++) {
        
        sx = 0.0; // line accumulator for max precision
        
        for(int j = 0; j< width; j++) {
            
            e = I0[i][j] - I1[i][j];
            sx += e * e;
        }
        s += sx;
    }
    return s;
}
// ---------------------------------------------------------------------------------
double of_calcMeanSquareError_F32(float32 **I0, float32 **I1, int height, int width)
// ---------------------------------------------------------------------------------
{
    // compute Mean Square Error (MSE)
    
    double se = of_calcSquareError_F32(I0, I1, height, width);
    double mse = se / (height * width);
    
    return mse;
}
// -----------------------------------------------------------------------------
double of_calcSquareError_U16(uint16 **I0, uint16 **I1, int height, int width)
// -----------------------------------------------------------------------------
{
    double e, sx, s; // pixel, line, total error
    
    s = 0.0;
    
    for(int i = 0; i < height; i++) {
        
        sx = 0.0; // line accumulator for max precision
        
        for(int j = 0; j< width; j++) {
            
            e = I0[i][j] - I1[i][j];
            sx += e * e;
        }
        s += sx;
    }
    return s;
}
// -------------------------------------------------------------------------------
double of_calcMeanSquareError_U16(uint16 **I0, uint16 **I1, int height, int width)
// -------------------------------------------------------------------------------
{
    // compute Mean Square Error (MSE)
    
    double se = of_calcSquareError_U16(I0, I1, height, width);
    double mse = se / (height * width);
    
    return mse;
}
// ----------------------------------------------------------------------------------------
double of_calcSquareError_U16_F32(uint16 **I0, int qi, float32 **I1, int height, int width)
// ----------------------------------------------------------------------------------------
{
    // attention qi code la difference de quantification entre I16 et F32
    double Qi = (double) (1 << qi);
    
    double s = 0.0;
    
    for(int i = 0; i < height; i++) {
        
        double sx = 0.0; // line accumulator for max precision
        
        for(int j = 0; j< width; j++) {
            double e = I0[i][j]/Qi - I1[i][j]; // pixel-difference
            sx += e * e; // line difference
        }
        s += sx; // total difference
    }
    return s;
}
// --------------------------------------------------------------------------------------------
double of_calcMeanSquareError_U16_F32(uint16 **I0, int qi, float32 **I1, int height, int width)
// --------------------------------------------------------------------------------------------
{
    // compute Mean Square Error (MSE)
    
    double se = of_calcSquareError_U16_F32(I0, qi, I1, height, width);
    double mse = se / (height * width);
    
    return mse;
}
// ----------------------------------------------------------------------------------------
double of_calcSquareError_S16_F32(sint16 **I0, int qi, float32 **I1, int height, int width)
// ----------------------------------------------------------------------------------------
{
    // attention qi code la difference de quantification entre I16 et F32
    double Qi = (double) (1 << qi);
    
    double s = 0.0;
    
    for(int i = 0; i < height; i++) {
        
        double sx = 0.0; // line accumulator for max precision
        
        for(int j = 0; j< width; j++) {
            double e = I0[i][j]/Qi - I1[i][j]; // pixel-difference
            sx += e * e; // line difference
        }
        s += sx; // total difference
    }
    return s;
}
// --------------------------------------------------------------------------------------------
double of_calcMeanSquareError_S16_F32(sint16 **I0, int qi, float32 **I1, int height, int width)
// --------------------------------------------------------------------------------------------
{
    // compute Mean Square Error (MSE)
    
    double se = of_calcSquareError_S16_F32(I0, qi, I1, height, width);
    double mse = se / (height * width);
    
    return mse;
}
// ----------------------------------------------------------------------------------------
double of_calcSquareError_S32_F32(sint32 **I0, int qi, float32 **I1, int height, int width)
// ----------------------------------------------------------------------------------------
{
    // attention qi code la difference de quantification entre I16 et F32
    double Qi = (double) (1 << qi);
    
    double s = 0.0;
    
    for(int i = 0; i < height; i++) {
        
        double sx = 0.0; // line accumulator for max precision
        
        for(int j = 0; j< width; j++) {
            double e = I0[i][j]/Qi - I1[i][j]; // pixel-difference
            sx += e * e; // line difference
        }
        s += sx; // total difference
    }
    return s;
}
// --------------------------------------------------------------------------------------------
double of_calcMeanSquareError_S32_F32(sint32 **I0, int qi, float32 **I1, int height, int width)
// --------------------------------------------------------------------------------------------
{
    // compute Mean Square Error (MSE)
    
    double se = of_calcSquareError_S32_F32(I0, qi, I1, height, width);
    double mse = se / (height * width);
    
    return mse;
}
// -------------------------------------------------------------------------------------
double of_calcRootMeanSquareError_F32(float32 **I0, float32 **I1, int height, int width)
// -------------------------------------------------------------------------------------
{
    // compute Root-Mean Square Error (SQUARE-ROOT SE)
    
    double mse  = of_calcMeanSquareError_F32(I0, I1, height, width);
    double rmse = sqrt(mse);
    
    return rmse;
}
// ------------------------------------------------------------------------------------
double of_calcInterpolationError_F32(float32 **I0, float32 **I1, int height, int width)
// ------------------------------------------------------------------------------------
{
    // IE = RMS error
    
    return of_calcRootMeanSquareError_F32(I0, I1, height, width);
}
// --------------------------------------------------------------------------------
double og_calcRelativeError_F32(float32 **Xref, float32 **X, int height, int width)
// --------------------------------------------------------------------------------
{
    double xref, x, s;
    
    s = 0.0;
    int c = 0;
    
    for(int i = 0; i < height; i++) {
        
        for(int j = 0; j< width; j++) {
            
            xref = Xref[i][j];
            x    = X[i][j];
          
            if(fabs(xref) > 0) {
                s += fabs(x-xref)/fabs(xref);
                c++;
            }
        }
    }
    
    if(c) {
        s /= c;
    } else {
        s = 0;
    }
    
    return s;
}