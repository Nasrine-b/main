/* ---------------- */
/* --- interp.c --- */
/* ---------------- */

/*
* interpolation functions
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
#include "interp.h"
#include "interp_operator.h"
// -------------------------------------------------------------------------------------------------
void interpNearest_F32(float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec)
// -------------------------------------------------------------------------------------------------
{
    #ifdef OPENMP
    #pragma omp parallel for firstprivate(height, width) shared(I, U, V, Irec)
    #endif
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {

            int iu, iv;
            float32 u, v;
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floorf(u);
            iv = (int) floorf(v);
                        
            Irec[i][j] = I[i+iv][j+iu];
        }
    }
}
// ------------------------------------------------------------------------------------------------
void interpLinear_F32(float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec)
// ------------------------------------------------------------------------------------------------
{
    #ifdef OPENMP
    #pragma omp parallel for firstprivate(height, width) shared(I, U, V, Irec)
    #endif
    
    for(int i = 0; i < height; i++) {
        for(int  j =0; j < width; j++) {
            
            int iu, iv;
            float32 u, v, fu, fv;
            
            float32 a00, a01;
            float32 a10, a11;
            
            float32 dx;
            float32 dy;
            
            float32 f, f0, f1;
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floorf(u);
            iv = (int) floorf(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            // x-computation
            dx = fu;
    
            // y-computation
            dy = fv;
            
            a00 = I[i+iv+0][j+iu+0]; a01 = I[i+iv+0][j+iu+1];
            a10 = I[i+iv+1][j+iu+0]; a11 = I[i+iv+1][j+iu+1];
            
            f0 = linear_interpolation(a00, a01, dx);
            f1 = linear_interpolation(a10, a11, dx);
            
            f0 = linear_interpolation(f0, f1, dy);
            
            Irec[i][j] = f;
        }
    }
}
#define CUBINTERP3(x1, x2, x3, x4) (   -x1 + 3 * x2 - 3 * x3 + x4)
#define CUBINTERP2(x1, x2, x3, x4) (2 * x1 - 5 * x2 + 4 * x3 - x4)
#define CUBINTERP1(x1, x2, x3, x4) (x1 - x3)
#define CUBINTERP0(x1, x2, x3, x4) x2

#define CUBINTERP(x1, x2, x3, x4, c3, c2, c1, c0) \
c3 = CUBINTERP3(x1, x2, x3, x4); \
c2 = CUBINTERP2(x1, x2, x3, x4); \
c1 = CUBINTERP1(x1, x2, x3, x4); \
c0 = CUBINTERP2(x1, x2, x3, x4)

#define HORNER3(c3, c2, c1, c0, x) (((c3 * x + c2) * x + c1) * x + c0)

// ---------------------------------------------------------------------------------------------------
void interpCubicPoly_F32(float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec)
// ---------------------------------------------------------------------------------------------------
{
    #ifdef OPENMP
    #pragma omp parallel for firstprivate(height, width) shared(I, U, V, Irec)
    #endif

    //puts("[interp]: debut");
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
        
            float32 u, v, fu, fv;
            
            float32 a11, a12, a13, a14;
            float32 a21, a22, a23, a24;
            float32 a31, a32, a33, a34;
            float32 a41, a42, a43, a44;
            
            float32 dx, dx2, dx3;
            float32 dy, dy2, dy3;
            
            float32 cx1, cx2, cx3, cx4;
            float32 cy1, cy2, cy3, cy4;
            
            float32 f1, f2, f3, f4;
            float32 f;
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            int iu = (int) floorf(u);
            int iv = (int) floorf(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            // x-computation
            dx  = fu;
            dx2 = dx * dx;
            dx3 = dx2 * dx;
            
            // y-computation
            dy  = fv;
            dy2 = dy * dy;
            dy3 = dy2 * dy;
            
            cx1 = -0.5f * (    dx3 - 2 * dx2 + dx);
            cx2 =  0.5f * (3 * dx3 - 5 * dx2 +  2);
            cx3 = -0.5f * (3 * dx3 - 4 * dx2 - dx);
            cx4 =  0.5f * (    dx3 -     dx2     );
            
            cy1 = -0.5f * (    dy3 - 2 * dy2 + dy);
            cy2 =  0.5f * (3 * dy3 - 5 * dy2 +  2);
            cy3 = -0.5f * (3 * dy3 - 4 * dy2 - dy);
            cy4 =  0.5f * (    dy3 -     dy2     );
            
            /*if( (i==1) && (j==0)) {
                
                printf("%10.4f %10.4f %3d %3d %10.4f %10.4f \n\n", u, v, iu, iv, fu, fv);
                
                printf("%10.4f %10.4f %10.4f %10.4f\n", cx1, cx2, cx3, cx4);
                printf("%10.4f %10.4f %10.4f %10.4f\n\n", cy1, cy2, cy3, cy4);
            }*/
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            f = cy1 * f1 + cy2 * f2 + cy3 * f3 + cy4 * f4;
            
            if( (i==1) && (j==0)) {
            
                /*printf("%10.4f %10.4f %10.4f %10.4f\n", a11, a12, a13, a14);
                printf("%10.4f %10.4f %10.4f %10.4f\n", a21, a22, a23, a24);
                printf("%10.4f %10.4f %10.4f %10.4f\n", a31, a32, a33, a34);
                printf("%10.4f %10.4f %10.4f %10.4f\n\n", a41, a42, a43, a44);
                
                printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n", f1, f2, f3, f4, f);*/
            }
            
            Irec[i][j] = f;
        }
    }
    //puts("[interp]: fin");
}
// ------------------------------------------------------------------------------------------------
void interpCubicU_F32(float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec)
// ------------------------------------------------------------------------------------------------
{
    #ifdef OPENMP
    #pragma omp parallel for firstprivate(height, width) shared(I, U, V, Irec)
    #endif

    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
        
            int iu, iv;
            float32 u, v, fu, fv;
            
            float32 a11, a12, a13, a14;
            float32 a21, a22, a23, a24;
            float32 a31, a32, a33, a34;
            float32 a41, a42, a43, a44;
            
            float32 dx;
            float32 dy;
            
            float32 c13,c12,c11,c10;
            float32 c23,c22,c21,c20;
            float32 c33,c32,c31,c30;
            float32 c43,c42,c41,c40;
            float32 c3, c2, c1, c0;
            
            float32 f1, f2, f3, f4;
            float32 f;
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floorf(u);
            iv = (int) floorf(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            // x-computation
            dx = fu;            
            // y-computation
            dy = fv;
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
                        
            // coeff d'interpolation cubique sans 1/2 et forme de Horner
            CUBINTERP(a11, a12, a13, a14, c13, c12, c11, c10); f1 = HORNER3(c13, c12, c11, c10, dx);
            CUBINTERP(a21, a22, a23, a24, c23, c22, c21, c20); f2 = HORNER3(c23, c22, c21, c20, dx);
            CUBINTERP(a31, a32, a33, a34, c33, c32, c31, c30); f3 = HORNER3(c33, c32, c31, c30, dx);
            CUBINTERP(a41, a42, a43, a44, c43, c42, c41, c40); f4 = HORNER3(c43, c42, c41, c40, dx);
            
            CUBINTERP(f1, f2, f3, f4, c3, c2, c1, c0); f = 0.25f * HORNER3(c3, c2, c1, c0, dy);
            
            Irec[i][j] = f;
        }
    }
}
// ---------------------------------------------------------------------------------------------------
void interpCubicCell_F32(float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec)
// ---------------------------------------------------------------------------------------------------
{
#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width) shared(I, U, V, Irec)
#endif
    for(int i = 0; i < height; i++) {
    
        // clipping en i
        float32 di0 = (float32) (0.0f   - i);
        float32 di1 = (float32) (height - i);
        
        for(int j = 0; j < width; j++) {
            
            // clipping en j
            float32 dj0 = (float32) (0.0f  - j);
            float32 dj1 = (float32) (width - j);
            
            int iu, iv;
            float32 u, v, fu, fv;
            
            float32 a11, a12, a13, a14;
            float32 a21, a22, a23, a24;
            float32 a31, a32, a33, a34;
            float32 a41, a42, a43, a44;
            
            float32 dx;
            float32 dy;
            
            float32 f1, f2, f3, f4;
            float32 f;
            
            u = U[i][j];
            v = V[i][j];
            
            // clipping
            /*if(u < (0.0f   - j)) u = (0.0f   - j);
            if(u > (width  - j)) u = (width  - j);
            if(v < (0.0f   - i)) v = (0.0f   - i);
            if(v > (height - i)) v = (height - i);*/
            
            // clipping
            if(u < dj0) u = dj0;
            if(u > dj1) u = dj1;
            if(v < di0) v = di0;
            if(v > di1) v = di1;
            
            // integer part
            iu = (int) floorf(u);
            iv = (int) floorf(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            // x-computation
            dx = fu;
            // y-computation
            dy = fv;
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
            
            f1 = cubic_interpolation_cell(a11, a12, a13, a14, dx);
            f2 = cubic_interpolation_cell(a21, a22, a23, a24, dx);
            f3 = cubic_interpolation_cell(a31, a32, a33, a34, dx);
            f4 = cubic_interpolation_cell(a41, a42, a43, a44, dx);
            
            f = cubic_interpolation_cell(f1, f2, f3, f4, dy);
            
            Irec[i][j] = f;
        }
    }
}
// -----------------------------------------------------------------------------------------------
void interpCubic_F32(float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec)
// -----------------------------------------------------------------------------------------------
{
    //interpCubicU_F32(I, U, V, height, width, Irec); // bug
    //interpCubicPoly_F32(I, U, V, height, width, Irec); // OK
    interpCubicCell_F32(I, U, V, height, width, Irec); // OK
}
// ------------------------------------------------------------------------------------------------------------
void gradientCubic_F32(float32 **I, float32 **U, float32 **V,int height, int width, float32 **Ix, float32 **Iy)
// ------------------------------------------------------------------------------------------------------------
{
    #ifdef OPENMP
    #pragma omp parallel for firstprivate(height, width) shared(I, U, V, Ix, Iy)
    #endif

    for(int i = 1; i < height; i++) {
        for(int j = 1; j < width; j++) {
        
            int iu, iv;
            float32 u, v, fu, fv;
            
            float32 a11, a12, a13, a14;
            float32 a21, a22, a23, a24;
            float32 a31, a32, a33, a34;
            float32 a41, a42, a43, a44;
            
            float32 dx, dx2, dx3;
            float32 dy, dy2, dy3;
            
            float32 cx1, cx2, cx3, cx4;
            float32 cy1, cy2, cy3, cy4;
            
            float32 gx1, gx2, gx3, gx4;
            float32 gy1, gy2, gy3, gy4;
            
            //float32 f;
            float32 f1, f2, f3, f4;
            float32 gx, gy;
            //float32 ggx, ggy; // pour debug
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floorf(u);
            iv = (int) floorf(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            // x-computation
            dx = fu;
            dx2 = dx * dx;
            dx3 = dx2 * dx;
            
            // y-computation
            dy = fv;
            dy2 = dy * dy;
            dy3 = dy2 * dy;
            
            // interpolation
            cx1 = -0.5f * (    dx3 - 2 * dx2 + dx);
            cx2 =  0.5f * (3 * dx3 - 5 * dx2 + 2 );
            cx3 = -0.5f * (3 * dx3 - 4 * dx2 - dx);
            cx4 =  0.5f * (    dx3 -     dx2     );
            
            cy1 = -0.5f * (    dy3 - 2 * dy2 + dy);
            cy2 =  0.5f * (3 * dy3 - 5 * dy2 + 2);
            cy3 = -0.5f * (3 * dy3 - 4 * dy2 - dy);
            cy4 =  0.5f * (    dy3 -     dy2);
            
            // gradient
            gx1 = -0.5f * (3 * dx2 -  4 * dx + 1);
            gx2 =  0.5f * (9 * dx2 - 10 * dx    );
            gx3 = -0.5f * (9 * dx2 -  8 * dx - 1);
            gx4 =  0.5f * (3 * dx2 -  2 * dx    );
            
            gy1 = -0.5f * (3 * dy2 -  4 * dy + 1);
            gy2 =  0.5f * (9 * dy2 - 10 * dy    );
            gy3 = -0.5f * (9 * dy2 -  8 * dy - 1);
            gy4 =  0.5f * (3 * dy2 -  2 * dy    );
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
            
            // formules 1D separables
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            gy = gy1 * f1 + gy2 * f2 + gy3 * f3 + gy4 * f4;
            
            f1 = gx1 * a11 + gx2 * a12 + gx3 * a13 + gx4 * a14;
            f2 = gx1 * a21 + gx2 * a22 + gx3 * a23 + gx4 * a24;
            f3 = gx1 * a31 + gx2 * a32 + gx3 * a33 + gx4 * a34;
            f4 = gx1 * a41 + gx2 * a42 + gx3 * a43 + gx4 * a44;
            
            gx = cy1 * f1 + cy2 * f2 + cy3 * f3 + cy4 * f4;
            
            Ix[i][j] = gx;
            Iy[i][j] = gy;
        }
    }
}
