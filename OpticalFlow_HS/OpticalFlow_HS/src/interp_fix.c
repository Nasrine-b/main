/* -------------------- */
/* --- interp_fix.c --- */
/* -------------------- */

/*
 * fonctions d'interpolations 2D
*/

/*
 * Copyright (c) 2006-2017 Lionel Lacassagne, all rights reserved
 * IEF, U-Psud, CNRS
 * LRI, U-Psud, CNRS
 * LIP6, UPMC, CNRS
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "nrc.h"
#include "nrdef.h"
#include "nrtype.h"

#include "of_macro.h"

#include "interp_fix.h"
#include "interp_operator.h"
#include "interp_operator_fix.h"

// ----------------------------------------------------------------------------------------------------------------
void interpLinear_I16_I32(uint16 **I, int qi, sint32 **U, sint32 **V, int qs, int height, int width, uint16 **Irec)
// ----------------------------------------------------------------------------------------------------------------
{
    sint64 Qi = CALC_Q(qi);
    sint64 Qs = CALC_Q(qs);
    
#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width, qi, qs, Qi, Qs) shared(I, U, V, Irec)
#endif
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            
            int iu, iv;
            int u, v, fu, fv;
            
            sint64 a00, a01;
            sint64 a10, a11;
            
            sint64 dx;
            sint64 dy;
            
            sint64 e0, e1;
            sint64 e;
            
            u = (int) U[i][j];
            v = (int) V[i][j];
            
            // integer part : (u) / Q et pas (u + Q/2) / Q car troncature
            iu = (int) u / Qs;
            iv = (int) v / Qs;
            
            // fraction part
            fu = u % Qs;
            fv = v % Qs;
            
            // x-computation
            dx = fu;
            // y-computation
            dy = fv;
            
            a00 = I[i+iv+0][j+iu+0]; a01 = I[i+iv+0][j+iu+1];
            a10 = I[i+iv+1][j+iu+0]; a11 = I[i+iv+1][j+iu+1];
            
            e0 = linear_interpolation_I64(a00, a01, dx, qs);
            e1 = linear_interpolation_I64(a10, a11, dx, qs);
            
            e = linear_interpolation_I64(e0, e1, dy, qs);
            
            // clipping
            if(e < 0) e = 0;
            if(e > Qi) e = Qi-1;
            
            Irec[i][j] = e;
        }
    }
}
// -------------------------------------------------------------------------------------------------------------------
void interpCubicIpol_I16_I32(uint16 **I, int qi, sint32 **U, sint32 **V, int qs, int height, int width, uint16 **Irec)
// -------------------------------------------------------------------------------------------------------------------
{
    sint32 Qi = CALC_Q(qi);
    sint32 Qs = CALC_Q(qs);
    
#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width, qi, qs, Qi, Qs) shared(I, U, V, Irec)
#endif
    
    for(int i = 0; i < height; i++) {
    
        // clipping en i
        sint32 di0 = (0      - i) * Qs;
        sint32 di1 = (height - i) * Qs;
        
        for(int j = 0; j < width; j++) {
            
            // clipping en j
            sint32 dj0 = (0     - j) * Qs;
            sint32 dj1 = (width - j) * Qs;
            
            int iu, iv;
            int u, v, fu, fv;
            
            sint64 a11, a12, a13, a14;
            sint64 a21, a22, a23, a24;
            sint64 a31, a32, a33, a34;
            sint64 a41, a42, a43, a44;
            
            sint64 dx;
            sint64 dy;
            
            sint64 f1, f2, f3, f4;
            sint64 f;
            
            u = (int) U[i][j];
            v = (int) V[i][j];
            
            // clipping
            if(u < dj0) u = dj0;
            if(u > dj1) u = dj1;
            if(v < di0) v = di0;
            if(v > di1) v = di1;
            
            // integer part : (u) / Q et pas (u + Q/2) / Q car troncature
            iu = (int) u / Qs;
            iv = (int) v / Qs;
            
            // fraction part
            fu = u % Qs;
            fv = v % Qs;
            
            // x-computation
            dx = fu;
            // y-computation
            dy = fv;
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
            
            f1 = cubic_interpolation_ipol_I64(a11, a12, a13, a14, dx, qs);
            f2 = cubic_interpolation_ipol_I64(a21, a22, a23, a24, dx, qs);
            f3 = cubic_interpolation_ipol_I64(a31, a32, a33, a34, dx, qs);
            f4 = cubic_interpolation_ipol_I64(a41, a42, a43, a44, dx, qs);
            
            f = cubic_interpolation_ipol_I64(f1, f2, f3, f4, dy, qs);
            
            // clipping
            if(f < 0)   f = 0;
            if(f >= Qi) f = Qi-1;
            
            Irec[i][j] = f;
        }
    }
}
// -----------------------------------------------------------------------------------------------------------------------
void interpCubicIpol_I16_F32(uint16 **I, int qi, sint32 **U, sint32 **V, int qs, int height, int width, uint16 **Irec)
// -------------------------------------------------------------------------------------------------------------------
{
    sint32 Qi = CALC_Q(qi);
    sint32 Qs = CALC_Q(qs);
    
#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width, qi, qs, Qi, Qs) shared(I, U, V, Irec)
#endif
    
    for(int i = 0; i < height; i++) {
        
        // clipping en i
        float32 di0 = (0      - i) * Qs;
        float32 di1 = (height - i) * Qs;
        
        for(int j = 0; j < width; j++) {
            
            // clipping en j
            float32 dj0 = (0     - j) * Qs;
            float32 dj1 = (width - j) * Qs;
            
            int iu, iv;
            int u, v, fu, fv;
            //float32 u_f, v_f;
            
            float32 a11, a12, a13, a14;
            float32 a21, a22, a23, a24;
            float32 a31, a32, a33, a34;
            float32 a41, a42, a43, a44;
            
            float32 dx;
            float32 dy;
            
            float32 f1, f2, f3, f4;
            float32 f;
            uint32  e;
            
            u = (int) U[i][j];
            v = (int) V[i][j];
            
            // clipping
            if(u < dj0) u = dj0;
            if(u > dj1) u = dj1;
            if(v < di0) v = di0;
            if(v > di1) v = di1;
            
            // integer part
            iu = (int) u / Qs;
            iv = (int) v / Qs;
            
            // fraction part
            fu = (float32) (u % Qs);
            fv = (float32) (v % Qs);
            
            // x-computation
            dx = fu;
            // y-computation
            dy = fv;
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
            
            f1 = cubic_interpolation_ipol_I16_F32(a11, a12, a13, a14, dx, qs);
            f2 = cubic_interpolation_ipol_I16_F32(a21, a22, a23, a24, dx, qs);
            f3 = cubic_interpolation_ipol_I16_F32(a31, a32, a33, a34, dx, qs);
            f4 = cubic_interpolation_ipol_I16_F32(a41, a42, a43, a44, dx, qs);
            
            f = cubic_interpolation_ipol_I16_F32(f1, f2, f3, f4, dy, qs);
            
            // clipping
            if(f < 0.0f) f = 0.0f;
            if(f >= Qi)  f = Qi - 1.0f;
            
            e = (uint16) lroundf(f);
            
            Irec[i][j] = (uint16) e;
        }
    }
}
// ===============================================================================================================
void interpCubic_I16_I32(uint16 **I, int qi, sint32 **U, sint32 **V, int qs, int height, int width, uint16 **Irec)
// ===============================================================================================================
{
    interpCubicIpol_I16_I32(I, qi, U, V, qs, height, width, Irec);
}
// ===============================================================================================================
void interpCubic_I16_F32(uint16 **I, int qi, sint32 **U, sint32 **V, int qs, int height, int width, uint16 **Irec)
// ===============================================================================================================
{
    // input and output in fixed-point
    // computation in floating-point
    
    interpCubicIpol_I16_F32(I, qi, U, V, qs, height, width, Irec);
}
