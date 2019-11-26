/* ---------------------- */
/* --- interp_extra.c --- */
/* ---------------------- */

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

#include "sequence.h" // pour tsSequence
//#include "of.h" // pour OF qui contient tsSequence

#include "interp.h"

/* ------------------------------------------------ */
void initInterpCubicLUT_F32(float32 **lut, int32 size)
/* ------------------------------------------------ */
{
    int i;
    float32 dx, dx2, dx3;
    float32 c1, c2, c3, c4;
    float32 step = 1.0f / (float32) size;
    
    for(i=-size; i<=size; i++) {
        
        dx  = step * (float32) i;
        dx2 = dx*dx;
        dx3 = dx2*dx;
        
        c1 = -0.5f * (  dx3-2*dx2+dx);
        c2 =  0.5f * (3*dx3-5*dx2+2);
        c3 = -0.5f * (3*dx3-4*dx2-dx);
        c4 =  0.5f * (  dx3-  dx2);     
        
        lut[i][1] = c1;
        lut[i][2] = c2;
        lut[i][3] = c3;
        lut[i][4] = c4;
    }
}
/* -------------------------------------------------- */
void initGradientCubicLUT_F32(float32 **lut, int32 size)
/* -------------------------------------------------- */
{
    int i;
    float32 dx, dx2, dx3;
    float32 c1, c2, c3, c4;
    float32 step = 1.0f / (float32) size;
    
    for(i=-size; i<=size; i++) {
        
        dx  = step * (float32) i;
        dx2 = dx*dx;
        dx3 = dx2*dx;
        
        c1 = -0.5f * (3*dx2-4*dx+1);
        c2 =  0.5f * (9*dx2-10*dx);
        c3 = -0.5f * (9*dx2-8*dx-1);
        c4 =  0.5f * (3*dx2-2*dx);
        
        lut[i][1] = c1;
        lut[i][2] = c2;
        lut[i][3] = c3;
        lut[i][4] = c4;
    }
}
/* -------------------------------------------------------------------------------------- */
void interpLinear_I8(uint8 **I, int height, int width, float32 **U, float32 **V, uint8 **Irec)
/* -------------------------------------------------------------------------------------- */
{
    int i,j;
    
    int iu, iv;
    float32 u, v, fu, fv;
    
    float32 a00, a01;
    float32 a10, a11;
    
    float32 cx1, cx2;
    float32 cy1, cy2;
    
    float32 f, f1, f2;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            cx1 = (1.0f-fu); cx2 = fu; 
            cy1 = (1.0f-fv); cy2 = fv; 
            
            a00 = I[i+iv+0][j+iu+0]; a01 = I[i+iv+0][j+iu+1];
            a10 = I[i+iv+1][j+iu+0]; a11 = I[i+iv+1][j+iu+1];
            
            f1 = cx1 * a00 + cx2 * a01;
            f2 = cx1 * a10 + cx2 * a11;
            
            f = cy1 * f1 + cy2 * f2;
            
            Irec[i][j] = (uint8) f;
        }
    }
}
/* ------------------------------------------------------------------------------------------- */
void interpLinear_I32(sint32 **I, int height, int width, float32 **U, float32 **V, sint32 **Irec)
/* ------------------------------------------------------------------------------------------- */
{
    int i,j;
    
    int iu, iv;
    float32 u, v, fu, fv;
    
    float32 a00, a01;
    float32 a10, a11;
    
    float32 cx1, cx2;
    float32 cy1, cy2;
    
    float32 f, f1, f2;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            cx1 = (1.0f-fu); cx2 = fu; 
            cy1 = (1.0f-fv); cy2 = fv; 
            
            a00 = (float32) I[i+iv+0][j+iu+0]; a01 = (float32) I[i+iv+0][j+iu+1];
            a10 = (float32) I[i+iv+1][j+iu+0]; a11 = (float32) I[i+iv+1][j+iu+1];
            
            f1 = cx1 * a00 + cx2 * a01;
            f2 = cx1 * a10 + cx2 * a11;
            
            f = cy1 * f1 + cy2 * f2;
            
            Irec[i][j] = (sint32)f;
        }
    }
}
/* ------------------------------------------------------------------------------------- */
void interpLacas_I8(uint8 **I, int height, int width, float32 **U, float32 **V, uint8 **Irec)
/* ------------------------------------------------------------------------------------- */
{
    int i,j;
    int iu, iv;
    float32 u, v, fu, fv;
    
    float32 dx, dx2, dx3;
    float32 dy, dy2, dy3;
    float32 cx1, cx2, cx3, cx4;
    float32 cy1, cy2, cy3, cy4;
    float32 f1, f2, f3, f4;
    float32 f;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
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
            
            cx1 = 0.25f * (dx3-3*dx2+2);
            cx2 = 0.25f * (dx3-3*dx2+4);
            cx3 = 0.25f * (dx3-6*dx2+9*dx);
            cx4 = 0.25f * (dx3-9*dx2+24*dx-16);
            
            cy1 = 0.25f * (dy3-3*dy2+2);
            cy2 = 0.25f * (dy3-3*dy2+4);
            cy3 = 0.25f * (dy3-6*dy2+9*dy);
            cy4 = 0.25f * (dy3-9*dy2+24*dy-16);  
            
            f1 = cx1 * I[i+iv-1][j+iu-1] + cx2 * I[i+iv-1][j+iu] + cx3 * I[i+iv-1][j+iu+1] + cx4 * I[i+iv-1][j+iu+2];
            f2 = cx1 * I[i+iv+0][j+iu-1] + cx2 * I[i+iv+0][j+iu] + cx3 * I[i+iv+0][j+iu+1] + cx4 * I[i+iv+0][j+iu+2];
            f3 = cx1 * I[i+iv+1][j+iu-1] + cx2 * I[i+iv+1][j+iu] + cx3 * I[i+iv+1][j+iu+1] + cx4 * I[i+iv+1][j+iu+2];
            f4 = cx1 * I[i+iv+2][j+iu-1] + cx2 * I[i+iv+2][j+iu] + cx3 * I[i+iv+2][j+iu+1] + cx4 * I[i+iv+2][j+iu+2];
            
            f = cy1*f1 + cy2*f2 + cy3*f3 + cy4*f4;
                        
            Irec[i][j] = (uint8) f;
        }
        //puts("");
    }
    //puts("");
}
/* -------------------------------------------------------------------------------------------- */
void interpLacas_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **Irec)
/* -------------------------------------------------------------------------------------------- */
{
    int i,j;
    int iu, iv;
    float32 u, v, fu, fv;
    
    float32 x, y;
    float32 dx, dx2, dx3;
    float32 dy, dy2, dy3;
    
    float32 cx1, cx2, cx3, cx4;
    float32 cy1, cy2, cy3, cy4;
    float32 f1, f2, f3, f4;
    float32 f;
    
    float32 a11, a12, a13, a14;
    float32 a21, a22, a23, a24;
    float32 a31, a32, a33, a34;
    float32 a41, a42, a43, a44;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
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
            
            /*cx1 = 0.25 * (dx3-3*dx2+2);
            cx2 = 0.25 * (dx3-3*dx2+4);
            cx3 = 0.25 * (dx3-6*dx2+9*dx);
            cx4 = 0.25 * (dx3-9*dx2+24*dx-16);
            
              cy1 = 0.25 * (dy3-3*dy2+2);
              cy2 = 0.25 * (dy3-3*dy2+4);
              cy3 = 0.25 * (dy3-6*dy2+9*dy);
            cy4 = 0.25 * (dy3-9*dy2+24*dy-16);*/  
            
            x = fu; y = fv;
            
            cx1 = (x-1)*(x-1)*(x+2)/4;
            cx2 = (x+1)*(x-2)*(x-2)/4;
            cx3 = x*(x-3)*(x-3)/4;
            cx4 = (x-1)*(x-4)*(x-4)/4;
            
            cy1 = (y-1)*(y-1)*(y+2)/4;
            cy2 = (y+1)*(y-2)*(y-2)/4;
            cy3 = y*(y-3)*(y-3)/4;
            cy4 = (y-1)*(y-4)*(y-4)/4;
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            f = cy1*f1 + cy2*f2 + cy3*f3 + cy4*f4;
                        
            Irec[i][j] = f;
        }
    }
}
/* ------------------------------------------------------------------------------------- */
void interpCubic_I8(uint8 **I, int height, int width, float32 **U, float32 **V, uint8 **Irec)
/* ------------------------------------------------------------------------------------- */
{
    int i,j;
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
    float32 f1, f2, f3, f4;
    float32 f;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
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
            
            cx1 = -0.5f * (  dx3-2*dx2+dx);
            cx2 =  0.5f * (3*dx3-5*dx2+2);
            cx3 = -0.5f * (3*dx3-4*dx2-dx);
            cx4 =  0.5f * (  dx3-  dx2);
            
            cy1 = -0.5f * (  dy3-2*dy2+dy);
            cy2 =  0.5f * (3*dy3-5*dy2+2);
            cy3 = -0.5f * (3*dy3-4*dy2-dy);
            cy4 =  0.5f * (  dy3-  dy2);
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            f = cy1*f1 + cy2*f2 + cy3*f3 + cy4*f4;
            
            Irec[i][j] = (uint8) f;
        }
    }
}
/* ------------------------------------------------------------------------------------------ */
void interpCubic_I32(sint32 **I, int height, int width, float32 **U, float32 **V, sint32 **Irec)
/* ------------------------------------------------------------------------------------------ */
{
    int i,j;
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
    float32 f1, f2, f3, f4;
    float32 f;
       
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
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
            
            cx1 = -0.5f * (  dx3-2*dx2+dx);
            cx2 =  0.5f * (3*dx3-5*dx2+2);
            cx3 = -0.5f * (3*dx3-4*dx2-dx);
            cx4 =  0.5f * (  dx3-  dx2);
            
            cy1 = -0.5f * (  dy3-2*dy2+dy);
            cy2 =  0.5f * (3*dy3-5*dy2+2);
            cy3 = -0.5f * (3*dy3-4*dy2-dy);
            cy4 =  0.5f * (  dy3-  dy2);
            
            a11 = (float32) I[i+iv-1][j+iu-1]; a12 = (float32) I[i+iv-1][j+iu]; a13 = (float32) I[i+iv-1][j+iu+1]; a14 = (float32) I[i+iv-1][j+iu+2];
            a21 = (float32) I[i+iv  ][j+iu-1]; a22 = (float32) I[i+iv  ][j+iu]; a23 = (float32) I[i+iv  ][j+iu+1]; a24 = (float32) I[i+iv  ][j+iu+2];
            a31 = (float32) I[i+iv+1][j+iu-1]; a32 = (float32) I[i+iv+1][j+iu]; a33 = (float32) I[i+iv+1][j+iu+1]; a34 = (float32) I[i+iv+1][j+iu+2];
            a41 = (float32) I[i+iv+2][j+iu-1]; a42 = (float32) I[i+iv+2][j+iu]; a43 = (float32) I[i+iv+2][j+iu+1]; a44 = (float32) I[i+iv+2][j+iu+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            f = cy1*f1 + cy2*f2 + cy3*f3 + cy4*f4;
            
            Irec[i][j] = (sint32) f;
        }
    }
}
/* ---------------------------------------------------------------------------------------------------- */
void interpCubic_F32_F64_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **Irec)
/* ---------------------------------------------------------------------------------------------------- */
{
    int i,j;
    int iu, iv;
    float32 u, v, fu, fv;
    
    float64 a11, a12, a13, a14;
    float64 a21, a22, a23, a24;
    float64 a31, a32, a33, a34;
    float64 a41, a42, a43, a44;
    
    float64 dx, dx2, dx3;
    float64 dy, dy2, dy3;
    float64 cx1, cx2, cx3, cx4;
    float64 cy1, cy2, cy3, cy4;
    float64 f1, f2, f3, f4;
    float64 f;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
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
            
            cx1 = -0.5 * (  dx3-2*dx2+dx);
            cx2 =  0.5 * (3*dx3-5*dx2+2);
            cx3 = -0.5 * (3*dx3-4*dx2-dx);
            cx4 =  0.5 * (  dx3-  dx2);
            
            cy1 = -0.5 * (  dy3-2*dy2+dy);
            cy2 =  0.5 * (3*dy3-5*dy2+2);
            cy3 = -0.5 * (3*dy3-4*dy2-dy);
            cy4 =  0.5 * (  dy3-  dy2);
                    
            a11 = (float64) I[i+iv-1][j+iu-1]; a12 = (float64) I[i+iv-1][j+iu]; a13 = (float64) I[i+iv-1][j+iu+1]; a14 = (float64) I[i+iv-1][j+iu+2];
            a21 = (float64) I[i+iv  ][j+iu-1]; a22 = (float64) I[i+iv  ][j+iu]; a23 = (float64) I[i+iv  ][j+iu+1]; a24 = (float64) I[i+iv  ][j+iu+2];
            a31 = (float64) I[i+iv+1][j+iu-1]; a32 = (float64) I[i+iv+1][j+iu]; a33 = (float64) I[i+iv+1][j+iu+1]; a34 = (float64) I[i+iv+1][j+iu+2];
            a41 = (float64) I[i+iv+2][j+iu-1]; a42 = (float64) I[i+iv+2][j+iu]; a43 = (float64) I[i+iv+2][j+iu+1]; a44 = (float64) I[i+iv+2][j+iu+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            f = cy1*f1 + cy2*f2 + cy3*f3 + cy4*f4;
            
            Irec[i][j] = (float32) f;
        }
    }
}
/* -------------------------------------------------------------------------------------------- */
void interpCubic_F64(float32 **I, int height, int width, float32 **U, float32 **V, float32 **Irec)
/* -------------------------------------------------------------------------------------------- */
{
    int i,j;
    int iu, iv;
    float64 u, v, fu, fv;
    
    float32 a11, a12, a13, a14;
    float32 a21, a22, a23, a24;
    float32 a31, a32, a33, a34;
    float32 a41, a42, a43, a44;
    
    float64 dx, dx2, dx3;
    float64 dy, dy2, dy3;
    float64 cx1, cx2, cx3, cx4;
    float64 cy1, cy2, cy3, cy4;
    float64 f1, f2, f3, f4;
    float64 f;
    
    #ifdef OPENMP
    #pragma omp parallel for private(i,j)
    #endif

    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
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
            
            cx1 = -0.5f * (  dx3-2*dx2+dx);
            cx2 =  0.5f * (3*dx3-5*dx2+2);
            cx3 = -0.5f * (3*dx3-4*dx2-dx);
            cx4 =  0.5f * (  dx3-  dx2);
            
            cy1 = -0.5f * (  dy3-2*dy2+dy);
            cy2 =  0.5f * (3*dy3-5*dy2+2);
            cy3 = -0.5f * (3*dy3-4*dy2-dy);
            cy4 =  0.5f * (  dy3-  dy2);
            
            a11 = (float32) I[i+iv-1][j+iu-1]; a12 = (float32) I[i+iv-1][j+iu]; a13 = (float32) I[i+iv-1][j+iu+1]; a14 = (float32) I[i+iv-1][j+iu+2];
            a21 = (float32) I[i+iv  ][j+iu-1]; a22 = (float32) I[i+iv  ][j+iu]; a23 = (float32) I[i+iv  ][j+iu+1]; a24 = (float32) I[i+iv  ][j+iu+2];
            a31 = (float32) I[i+iv+1][j+iu-1]; a32 = (float32) I[i+iv+1][j+iu]; a33 = (float32) I[i+iv+1][j+iu+1]; a34 = (float32) I[i+iv+1][j+iu+2];
            a41 = (float32) I[i+iv+2][j+iu-1]; a42 = (float32) I[i+iv+2][j+iu]; a43 = (float32) I[i+iv+2][j+iu+1]; a44 = (float32) I[i+iv+2][j+iu+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            f = cy1*f1 + cy2*f2 + cy3*f3 + cy4*f4;
            
            Irec[i][j] = (float32) f;
        }
    }
}
/* ---------------------------------------------------------------------------------------------------------------------------- */
void interpCubicLUT_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **lutI, int32 sizeI, float32 **Irec)
/* ---------------------------------------------------------------------------------------------------------------------------- */
{
    int i,j;
    int iu, iv;
    float32 u, v, fu, fv;
    
    float32 a11, a12, a13, a14;
    float32 a21, a22, a23, a24;
    float32 a31, a32, a33, a34;
    float32 a41, a42, a43, a44;
    
    float32 dx;
    float32 dy;
    float32 szI = (float32) sizeI;
    int     indx, indy;
    
    float32 cx1, cx2, cx3, cx4;
    float32 cy1, cy2, cy3, cy4;
    float32 f1, f2, f3, f4;
    float32 f;

    #ifdef OPENMP
    #pragma omp parallel for private(i,j)
    #endif

    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            // LUT addresses
            dx = fu;
            dy = fv;
            
            indx = (int) (dx * szI);
            indy = (int) (dy * szI);
            
            // coeff extraction
            
            cx1 = lutI[indx][1];
            cx2 = lutI[indx][2];
            cx3 = lutI[indx][3];
            cx4 = lutI[indx][4];
            
            cy1 = lutI[indy][1];
            cy2 = lutI[indy][2];
            cy3 = lutI[indy][3];
            cy4 = lutI[indy][4];
            
            a11 = (float32) I[i+iv-1][j+iu-1]; a12 = (float32) I[i+iv-1][j+iu]; a13 = (float32) I[i+iv-1][j+iu+1]; a14 = (float32) I[i+iv-1][j+iu+2];
            a21 = (float32) I[i+iv  ][j+iu-1]; a22 = (float32) I[i+iv  ][j+iu]; a23 = (float32) I[i+iv  ][j+iu+1]; a24 = (float32) I[i+iv  ][j+iu+2];
            a31 = (float32) I[i+iv+1][j+iu-1]; a32 = (float32) I[i+iv+1][j+iu]; a33 = (float32) I[i+iv+1][j+iu+1]; a34 = (float32) I[i+iv+1][j+iu+2];
            a41 = (float32) I[i+iv+2][j+iu-1]; a42 = (float32) I[i+iv+2][j+iu]; a43 = (float32) I[i+iv+2][j+iu+1]; a44 = (float32) I[i+iv+2][j+iu+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            f = cy1*f1 + cy2*f2 + cy3*f3 + cy4*f4;
            
            Irec[i][j] = f;
        }
    }
}
/* ----------------------------------------------------------------------- */
void sobel_I8F32(uint8 **I, int height, int width, float32 **Ix, float32 **Iy)
/* ----------------------------------------------------------------------- */
{
    int i,j;
    
    sint16 a0,a1,a2;
    sint16 b0,b1,b2;
    sint16 c0,c1,c2;
    
    sint32 gx, gy;
    
    for(i=0; i<height; i++) {
        for(j=0; j<=width; j++) {
            
            a0 = I[i-1][j-1]; a1 = I[i-1][j]; a2 = I[i-1][j+1];
            b0 = I[i  ][j-1]; b1 = I[i  ][j]; b2 = I[i  ][j+1];
            c0 = I[i+1][j-1]; c1 = I[i+1][j]; c2 = I[i+1][j+1];
            
            gx = (a2-a0) + 2*(b2-b0) + (c2-c0);
            gy = (c0-a0) + 2*(c1-a1) + (c2-a2);
            
            //gx = gx/4;
            //gy = gy/4;
            
            Ix[i][j] = (float32) gx;
            Iy[i][j] = (float32) gy;
        }
    }
}
/* ------------------------------------------------------------------------ */
void sobel_F32(float32 **I, int height, int width, float32 **Ix, float32 **Iy)
/* ------------------------------------------------------------------------ */
{
    int i,j;
    
    float32 a0,a1,a2;
    float32 b0,b1,b2;
    float32 c0,c1,c2;
    
    float32 gx, gy;
    
    for(i=0; i<height; i++) {
        for(j=0; j<=width; j++) {
            
            a0 = I[i-1][j-1]; a1 = I[i-1][j]; a2 = I[i-1][j+1];
            b0 = I[i  ][j-1]; b1 = I[i  ][j]; b2 = I[i  ][j+1];
            c0 = I[i+1][j-1]; c1 = I[i+1][j]; c2 = I[i+1][j+1];
            
            gx = (a2-a0) + 2*(b2-b0) + (c2-c0);
            gy = (c0-a0) + 2*(c1-a1) + (c2-a2);
            
            //gx = gx/4;
            //gy = gy/4;
            
            Ix[i][j] = gx;
            Iy[i][j] = gy;
        }
    }
}
/* -------------------------------------------------------------------------- */
void gradient_I8F32(uint8 **I, int height, int width, float32 **Ix, float32 **Iy)
/* -------------------------------------------------------------------------- */
{
    sobel_I8F32(I, height, width, Ix, Iy);
}
/* ------------------------------------------------------------------------------------------------------- */
void gradientLinear_I8(uint8 **I, int height, int width, float32 **U, float32 **V, float32 **Ix, float32 **Iy)
/* ------------------------------------------------------------------------------------------------------- */
{
    int i,j;
    
    int iu, iv;
    float32 u, v, fu, fv;
    
    float32 a00, a01;
    float32 a10, a11;
    
    //float32 cx1, cx2;
    //float32 cy1, cy2;
    
    float32 gx, gx1, gx2;
    float32 gy, gy1, gy2;
    
    //float32 f, f1, f2;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            gx1 = (1.0f-fv); gx2 = fv; // attention u et v sont inverses
            gy1 = (1.0f-fu); gy2 = fu; 
            
            a00 = I[i+iv+0][j+iu+0]; a01 = I[i+iv+0][j+iu+1];
            a10 = I[i+iv+1][j+iu+0]; a11 = I[i+iv+1][j+iu+1];
            
            gx = gx1 * (a01-a00) + gx2 * (a11-a10);
            gy = gy1 * (a10-a00) + gy2 * (a11-a01);
            
            Ix[i][j] = (uint8) gx;
            Iy[i][j] = (uint8) gy;
        }
    }
}
/* ----------------------------------------------------------------------------------------------------------- */
void gradientLinear_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **Ix, float32 **Iy)
/* ----------------------------------------------------------------------------------------------------------- */
{
    int i,j;
    
    int iu, iv;
    float32 u, v, fu, fv;
    
    float32 a00, a01;
    float32 a10, a11;
    
    //float32 cx1, cx2;
    //float32 cy1, cy2;
    
    float32 gx, gx1, gx2;
    float32 gy, gy1, gy2;
    
    //float32 f, f1, f2;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            gx1 = (1.0f-fv); gx2 = fv; // attention u et v sont inverses
            gy1 = (1.0f-fu); gy2 = fu; 
            
            a00 = I[i+iv+0][j+iu+0]; a01 = I[i+iv+0][j+iu+1];
            a10 = I[i+iv+1][j+iu+0]; a11 = I[i+iv+1][j+iu+1];
            
            gx = gx1 * (a01-a00) + gx2 * (a11-a10);
            gy = gy1 * (a10-a00) + gy2 * (a11-a01);
            
            Ix[i][j] = gx;
            Iy[i][j] = gy;
        }
    }
}
/* ------------------------------------------------------------------------------------------------------ */
void gradientCubic_I8(uint8 **I, int height, int width, float32 **U, float32 **V, float32 **Ix, float32 **Iy)
/* ------------------------------------------------------------------------------------------------------ */
{
    int i,j;
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
    
    float32 f1, f2, f3, f4;
    float32 gx, gy;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
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
            cx1 = -0.5f * (dx3-2*dx2+dx);
            cx2 =  0.5f * (3*dx3-5*dx2+2);
            cx3 = -0.5f * (3*dx3-4*dx2-dx);
            cx4 =  0.5f * (dx3-dx2);
            
            cy1 = -0.5f * (dy3-2*dy2+dy);
            cy2 =  0.5f * (3*dy3-5*dy2+2);
            cy3 = -0.5f * (3*dy3-4*dy2-dy);
            cy4 =  0.5f * (dy3-dy2);
            
            // gradient
            gx1 = -0.5f * (3*dx2-4*dx+1);
            gx2 =  0.5f * (9*dx2-10*dx);
            gx3 = -0.5f * (9*dx2-8*dx-1);
            gx4 =  0.5f * (3*dx2-2*dx);
            
            gy1 = -0.5f * (3*dy2-4*dy+1);
            gy2 =  0.5f * (9*dy2-10*dy);
            gy3 = -0.5f * (9*dy2-8*dy-1);
            gy4 =  0.5f * (3*dy2-2*dy);
            
            a11 = I[i+iv-1][j+iu-1]; a12 = I[i+iv-1][j+iu]; a13 = I[i+iv-1][j+iu+1]; a14 = I[i+iv-1][j+iu+2];
            a21 = I[i+iv  ][j+iu-1]; a22 = I[i+iv  ][j+iu]; a23 = I[i+iv  ][j+iu+1]; a24 = I[i+iv  ][j+iu+2];
            a31 = I[i+iv+1][j+iu-1]; a32 = I[i+iv+1][j+iu]; a33 = I[i+iv+1][j+iu+1]; a34 = I[i+iv+1][j+iu+2];
            a41 = I[i+iv+2][j+iu-1]; a42 = I[i+iv+2][j+iu]; a43 = I[i+iv+2][j+iu+1]; a44 = I[i+iv+2][j+iu+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            gy = gy1 * f1 + gy2 * f2 + gy3 * f3 * gy4 * f4;
            
            f1 = gx1 * a11 + gx2 * a12 + gx3 * a13 + gx4 * a14;
            f2 = gx1 * a21 + gx2 * a22 + gx3 * a23 + gx4 * a24;
            f3 = gx1 * a31 + gx2 * a32 + gx3 * a33 + gx4 * a34;
            f4 = gx1 * a41 + gx2 * a42 + gx3 * a43 + gx4 * a44;
            
            gx = cy1 * f1 + cy2 * f2 + cy3 * f3 * cy4 * f4;
            
            Ix[i][j] = gx;
            Iy[i][j] = gy;
        }
    }
}
/* --------------------------------------------------------------------------------- */
void gradientCubic0_F32(float32 **I, int height, int width, float32 **Ix, float32 **Iy)
/* --------------------------------------------------------------------------------- */
{
    int i,j;
    
    float32 a12, a21, a23, a32;
    
    float32 gx, gy;
        
    #ifdef OPENMP
    #pragma omp parallel for private(i,j,a12,a21,a23,a32,gx,gy)
    #endif

    for(i=1; i<height; i++) {
        for(j=1; j<width; j++) { 
                       
            a12 = I[i-1][j  ];
            a21 = I[i  ][j-1];
            a23 = I[i  ][j+1];
            a32 = I[i+1][j  ];
            
            gx = 0.5f * (a23-a21);
            gy = 0.5f * (a32-a12);
            
            Ix[i][j] = gx;
            Iy[i][j] = gy;
        }
    }
    return;
}
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
void gradientCubicLUT_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **lutI, int32 sizeI, float32 **lutG, int32 sizeG, float32 **Ix, float32 **Iy)
/* ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
{
    int i,j;
    int iu, iv;
    float32 u, v, fu, fv;
    
    float32 a11, a12, a13, a14;
    float32 a21, a22, a23, a24;
    float32 a31, a32, a33, a34;
    float32 a41, a42, a43, a44;
    
    float32 dx;
    float32 dy;
    float32 szI = (float32) sizeI;
    float32 szG = (float32) sizeG;
    int     indxI, indyI;
    int     indxG, indyG;
    
    float32 cx1, cx2, cx3, cx4;
    float32 cy1, cy2, cy3, cy4;
    
    float32 gx1, gx2, gx3, gx4;
    float32 gy1, gy2, gy3, gy4;
    
    //float32 f;
    float32 f1, f2, f3, f4;
    float32 gx, gy;
    //float32 ggx, ggy; // pour debug
    

    //#pragma omp parallel for private(i,j)
    for(i=1; i<height; i++) {
        for(j=1; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
            // fraction part
            fu = u - iu;
            fv = v - iv;
            
            // LUT addresses
            dx = fu;
            dy = fv;
            
            indxI = (int) (dx * szI);
            indyI = (int) (dy * szI);
            
            indxG = (int) (dx * szG);
            indyG = (int) (dy * szG);
            
            // interpolation
            cx1 = lutI[indxI][1];
            cx2 = lutI[indxI][2];
            cx3 = lutI[indxI][3];
            cx4 = lutI[indxI][4];
            
            cy1 = lutI[indyI][1];
            cy2 = lutI[indyI][2];
            cy3 = lutI[indyI][3];
            cy4 = lutI[indyI][4];
            
            // gradient
            gx1 = lutG[indxG][1];
            gx2 = lutG[indxG][2];
            gx3 = lutG[indxG][3];
            gx4 = lutG[indxG][4];
            
            gy1 = lutG[indyG][1];
            gy2 = lutG[indyG][2];
            gy3 = lutG[indyG][3];
            gy4 = lutG[indyG][4];
            
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
    return;
}
/* ------------------------------------------------------------------------------------------------------ */
void gradientLacas_I8(uint8 **I, int height, int width, float32 **U, float32 **V, float32 **Ix, float32 **Iy)
/* ------------------------------------------------------------------------------------------------------ */
{
    int i,j;
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
    
    float32 f1, f2, f3, f4;
    float32 gx, gy;
    
    for(i=0; i<height; i++) {
        for(j=0; j<width; j++) { 
            
            u = U[i][j];
            v = V[i][j];
            
            // integer part
            iu = (int) floor(u);
            iv = (int) floor(v);
            
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
            cx1 = 0.25f * (dx3-3*dx2+2);
            cx2 = 0.25f * (dx3-3*dx2+4);
            cx3 = 0.25f * (dx3-6*dx2+9*dx);
            cx4 = 0.25f * (dx3-9*dx2+24*dx-16);
            
            cy1 = 0.25f * (dy3-3*dy2+2);
            cy2 = 0.25f * (dy3-3*dy2+4);
            cy3 = 0.25f * (dy3-6*dy2+9*dy);
            cy4 = 0.25f * (dy3-9*dy2+24*dy-16);  
            
            // gradient
            gx1 = 0.75f * (dx2-1);
            gx2 = 0.75f * (dx2-dx);
            gx3 = 0.75f * (dx2-4*dx+3);
            gx4 = 0.75f * (dx2-6*dx+8);
            
            gy1 = 0.75f * (dy2-1);
            gy2 = 0.75f * (dy2-dy);
            gy3 = 0.75f * (dy2-4*dy+3);
            gy4 = 0.75f * (dy2-6*dy+8);
            
            a11 = I[i+iu-1][j+iv-1]; a12 = I[i+iu-1][j+iv]; a13 = I[i+iu-1][j+iv+1]; a14 = I[i+iu-1][j+iv+2];
            a21 = I[i+iu  ][j+iv-1]; a22 = I[i+iu  ][j+iv]; a23 = I[i+iu  ][j+iv+1]; a24 = I[i+iu  ][j+iv+2];
            a31 = I[i+iu+1][j+iv-1]; a32 = I[i+iu+1][j+iv]; a33 = I[i+iu+1][j+iv+1]; a34 = I[i+iu+1][j+iv+2];
            a41 = I[i+iu+2][j+iv-1]; a42 = I[i+iu+2][j+iv]; a43 = I[i+iu+2][j+iv+1]; a44 = I[i+iu+2][j+iv+2];
            
            f1 = cx1 * a11 + cx2 * a12 + cx3 * a13 + cx4 * a14;
            f2 = cx1 * a21 + cx2 * a22 + cx3 * a23 + cx4 * a24;
            f3 = cx1 * a31 + cx2 * a32 + cx3 * a33 + cx4 * a34;
            f4 = cx1 * a41 + cx2 * a42 + cx3 * a43 + cx4 * a44;
            
            gy = gy1 * f1 + gy2 * f2 + gy3 * f3 * gy4 * f4;
            
            f1 = gx1 * a11 + gx2 * a12 + gx3 * a13 + gx4 * a14;
            f2 = gx1 * a21 + gx2 * a22 + gx3 * a23 + gx4 * a24;
            f3 = gx1 * a31 + gx2 * a32 + gx3 * a33 + gx4 * a34;
            f4 = gx1 * a41 + gx2 * a42 + gx3 * a43 + gx4 * a44;
            
            gx = cy1 * f1 + cy2 * f2 + cy3 * f3 * cy4 * f4;
            
            Ix[i][j] = gx;
            Iy[i][j] = gy;
        }
    }
}
