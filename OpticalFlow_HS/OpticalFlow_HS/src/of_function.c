/* --------------------- */
/* --- OF_function.c --- */
/* --------------------- */

/*
 * fonctions generalistes d'OF
 */

/*
 * Copyright (c) 2006-2009 Lionel Lacassagne, all rights reserved, IEF, U-PSud, CNRS
 * Copyright (c) 2018-2018 Lionel Lacassagne, all rights reserved, LIP6, SU, CNRS
 */

/*
 * 2006-2009 creation
 * 2018: mise a jour (Burt)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>


#include "nrc.h"
#include "nrtype.h"
#include "nrdef.h"
#include "nralloc.h"
#include "nrarith.h"
#include "nrio.h"

#include "sequence.h"

#include "of_config.h"

#include "of_function.h"

#ifdef ENABLE_X
#include "of_average_fast.h"
#endif


// ----------------------
int is_nan_F32(float32 f)
// ----------------------
{
    if(f != f)
        return 1;
    else
        return 0;
}
// ----------------------
int is_nan_F64(float64 f)
// ----------------------
{
    if(f != f)
        return 1;
    else
        return 0;
}
// -------------------------------------------------------------------------
int isnan_f32matrix(float32 **m, int i0, int i1, int j0, int j1, char *name)
// -------------------------------------------------------------------------
{
    /*int i, j;
    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            if(is_nan_F32(m[i][j])) {
                printf("%s %d %d is a NaN\n", name, i, j);
            }
            if(isinf(m[i][j])) {
                printf("%s %d %d is inf\n", name, i, j);
            }
        }
    }
    puts("");*/
    return 0;
}
/* ------------------------------------------------------------------ */
void makeBorder_I8(uint8 **T, int i0, int i1, int j0, int j1, int border)
/* ------------------------------------------------------------------ */
{
    int i, j, k;

    // --------------- //
    // -- bord haut -- //
    // --------------- //

    for(k=1; k<=border; k++) {
        for(j=j0; j<=j1; j++) {
            T[i0-k][j] = T[i0][j];
        }
    }

    // -------------- //
    // -- bord bas -- //
    // -------------- //

    for(k=1; k<=border; k++) {
        for(j=j0; j<=j1; j++) {
            T[i1+k][j] = T[i1][j];
        }
    }

    // --------------------------- //
    // -- bords gauche et droit -- //
    // --------------------------- //

    for(i=i0-border; i<=i1+border; i++) {
        for(k=1; k<=border; k++) {
            T[i][j0-k] = T[i][j0];
            T[i][j1+k] = T[i][j1];
        }
    }
}
/* ------------------------------------------------------------------ */
void zeroBorder_I8(uint8 **T, int i0, int i1, int j0, int j1, int border)
/* ------------------------------------------------------------------ */
{
    int i, j, k;
    
    // --------------- //
    // -- bord haut -- //
    // --------------- //
    
    for(k=1; k<=border; k++) {
        for(j=j0; j<=j1; j++) {
            T[i0-k][j] = 0;
        }
    }
    
    // -------------- //
    // -- bord bas -- //
    // -------------- //
    
    for(k=1; k<=border; k++) {
        for(j=j0; j<=j1; j++) {
            T[i1+k][j] = 0;
        }
    }
    
    // --------------------------- //
    // -- bords gauche et droit -- //
    // --------------------------- //
    
    for(i=i0-border; i<=i1+border; i++) {
        for(k=1; k<=border; k++) {
            T[i][j0-k] = 0;
            T[i][j1+k] = 0;
        }
    }   
}
/* -------------------------------------- */
void extendDataZ2_I8(uint8 **X, int h, int w)
/* -------------------------------------- */
{
    zeroBorder_I8(X, 0, h-1, 0, w-1, 2);
}
/* -------------------------------------- */
void extendDataD2_I8(uint8 **X, int h, int w)
/* -------------------------------------- */
{
    makeBorder_I8(X, 0, h-1, 0, w-1, 2);
}
/* ---------------------------------------------------------------------- */
void makeBorder_F32(float32 **T, int i0, int i1, int j0, int j1, int border)
/* ---------------------------------------------------------------------- */
{
    int i, j, k;

    // --------------- //
    // -- bord haut -- //
    // --------------- //

    for(k=1; k<=border; k++) {
        for(j=j0; j<=j1; j++) {
            T[i0-k][j] = T[i0][j];
        }
    }

    // -------------- //
    // -- bord bas -- //
    // -------------- //

    for(k=1; k<=border; k++) {
        for(j=j0; j<=j1; j++) {
            T[i1+k][j] = T[i1][j];
        }
    }

    // --------------------------- //
    // -- bords gauche et droit -- //
    // --------------------------- //

    for(i=i0-border; i<=i1+border; i++) {
        for(k=1; k<=border; k++) {
            T[i][j0-k] = T[i][j0];
            T[i][j1+k] = T[i][j1];
        }
    }
}
/* ---------------------------------------------------------------------- */
void zeroBorder_F32(float32 **T, int i0, int i1, int j0, int j1, int border)
/* ---------------------------------------------------------------------- */
{
    int i, j, k;
    
    // --------------- //
    // -- bord haut -- //
    // --------------- //
    
    for(k=1; k<=border; k++) {
        for(j=j0; j<=j1; j++) {
            T[i0-k][j] = 0;
        }
    }
    
    // -------------- //
    // -- bord bas -- //
    // -------------- //
    
    for(k=1; k<=border; k++) {
        for(j=j0; j<=j1; j++) {
            T[i1+k][j] = 0;
        }
    }
    
    // --------------------------- //
    // -- bords gauche et droit -- //
    // --------------------------- //
    
    for(i=i0-border; i<=i1+border; i++) {
        for(k=1; k<=border; k++) {
            T[i][j0-k] = 0;
            T[i][j1+k] = 0;
        }
    }   
}
/* ------------------------------------------ */
void extendDataZ2_F32(float32 **X, int h, int w)
/* ------------------------------------------ */
{
    zeroBorder_F32(X, 0, h-1, 0, w-1, 2);
}
/* ------------------------------------------ */
void extendDataD2_F32(float32 **X, int h, int w)
/* ------------------------------------------ */
{
    // ERREUR DE CONCEPTION
    // R = 2
    // sauf si pour F uniquement (ie pour interpolation bi-cubique necessitant un rayon de 2)
    makeBorder_F32(X, 0, h-1, 0, w-1, 2);
}
/* ----------------------------------------------------------------------------------- */
float32 max_matrix_pos(float32 **X, int i0, int i1, int j0, int j1, int *ipos, int *jpos)
/* ----------------------------------------------------------------------------------- */
{
    int i, j;
    int ii, jj;
    float32 m;

    ii = i0;
    jj = j0;
    m = X[ii][jj];

    for(i=i0; i<=i1; i++) {
        for(j=j0; j<=j1; j++) {
            if(X[i][j] > m) {
                m = X[i][j];
                ii = i;
                jj = j;
            }
        }
    }
    *ipos = ii;
    *jpos = jj;
    return m;
}
// -----------------------------------------------------------------------------------------
void convert_f32matrix_ui8matrix_sat(float32 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------------
{
    float32 x;
    uint8 y;
    
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            
            x = (sint16) X[i][j];
            
            if     (x < 0  ) y = 0;
            else if(x > 255) y = 255;
            else             y = (uint8) x;
            
            Y[i][j] = y;
        }
    }
}
// -----------------------------------------------------------------------------------------
void initKernel_F32(float32 **K, float32 c0, float32 c1, float32 c2, float32 c3, float32 c4)
// -----------------------------------------------------------------------------------------
{
    //noyau separable
    
    K[-2][-2] = c0 * c0; K[-2][-1] = c0 * c1; K[-2][0] = c0 * c2; K[-2][1] = c0 * c3; K[-2][2] = c0 * c4;
    K[-1][-2] = c1 * c0; K[-1][-1] = c1 * c1; K[-1][0] = c1 * c2; K[-1][1] = c1 * c3; K[-1][2] = c1 * c4;
    K[ 0][-2] = c2 * c0; K[ 0][-1] = c2 * c1; K[ 0][0] = c2 * c2; K[ 0][1] = c2 * c3; K[ 0][2] = c2 * c4;
    K[ 1][-2] = c3 * c0; K[ 1][-1] = c3 * c1; K[ 1][0] = c3 * c2; K[ 1][1] = c3 * c3; K[ 1][2] = c3 * c4;
    K[ 2][-2] = c4 * c0; K[ 2][-1] = c4 * c1; K[ 2][0] = c4 * c2; K[ 2][1] = c4 * c3; K[ 2][2] = c4 * c4;
}
// --------------------------------------------
void initBurtKernel_F32(float32 a, float32 **K)
// --------------------------------------------
{
    float32 c0, c1, c2, c3, c4;

    c0 = (float) (1.0/4.0-a/2.0);
    c1 = (float) (1.0/4.0);
    c2 = a;
    c3 = (float) (1.0/4.0);
    c4 = (float) (1.0/4.0-a/2.0);

    initKernel_F32(K, c0, c1, c2, c3, c4);

}
/* ---------------------------- */
void initGauss5Kernel_F32(float32 **K)
/* ---------------------------- */
{
    float32 c0, c1, c2, c3, c4;

    c0 = (float) (1.0/16);
    c1 = (float) (4.0/16);
    c2 = (float) (6.0/16);
    c3 = (float) (4.0/16);
    c4 = (float) (1.0/16);

    initKernel_F32(K, c0, c1, c2, c3, c4);

}
/* -------------------------------- */
void initGauss3Kernel_F32(float32 **K)
/* -------------------------------- */
{
    float32 c0, c1, c2, c3, c4;

    c0 = (float) 0;
    c1 = (float) (1.0/4.0);
    c2 = (float) (2.0/4.0);
    c3 = (float) (1.0/4.0);
    c4 = (float) 0;

    initKernel_F32(K, c0, c1, c2, c3, c4);
}
/* ---------------------------------- */
void initAverage5Kernel_F32(float32 **K)
/* ---------------------------------- */
{
    float32 c0, c1, c2, c3, c4;

    c0 = (float)(1.0/5.0);
    c1 = (float)(1.0/5.0);
    c2 = (float)(1.0/5.0);
    c3 = (float)(1.0/5.0);
    c4 = (float)(1.0/5.0);

    initKernel_F32(K, c0, c1, c2, c3, c4);
}
/* ---------------------------------- */
void initAverage3Kernel_F32(float32 **K)
/* ---------------------------------- */
{
    float32 c0, c1, c2, c3, c4;

    c0 = (float) 0.0f;
    c1 = (float)(1.0f/3.0f);
    c2 = (float)(1.0f/3.0f);
    c3 = (float)(1.0f/3.0f);
    c4 = (float) 0.0f;

    initKernel_F32(K, c0, c1, c2, c3, c4);
}
/* -------------------------------------------------------------------------------- */
void convolution5x5_F32(float32 **X, int height, int width, float32 **K,  float32 **Y)
/* -------------------------------------------------------------------------------- */
{
    #ifdef OPENMP
    #pragma omp parallel for
    #endif

    for(int i=0; i<height; i++) {
        for(int j=0; j<width; j++) {

            float32 f = 0.0f;
            
            for(int k=-2; k<=+2; k++) {
                for(int l=-2; l<=+2; l++) {
                    f += X[i+k][j+l] * K[k][l];
                }
            }
            Y[i][j] = (uint8) f;
        }
    }
}
// -----------------------------------------------------------------------------------------------
void convolution5x5_decimation_I8_F32_I8(uint8 **X, int height, int width, float32 **K, uint8 **Y)
// -----------------------------------------------------------------------------------------------
{
#ifdef OPENMP
#pragma omp parallel for
#endif

    for(int i=0; i<height; i+=2) {
        for(int j=0; j<width; j+=2) {

            float32 f = 0.0f;
            
            for(int k=-2; k<=+2; k++) {
                for(int l=-2; l<=+2; l++) {
                    f += X[i+k][j+l] * K[k][l];
                }
            }
            
            Y[i/2][j/2] = (uint8) f;
        }
    }
}
/* --------------------------------------------------------------------------------------------------- */
void convolution5x5_decimation_F32_F32_F32(float32 **X, int height, int width, float32 **K,  float32 **Y)
/* --------------------------------------------------------------------------------------------------- */
{
    #ifdef OPENMP
    #pragma omp parallel for
    #endif

    for(int i=0; i<height; i+=2) {
        for(int j=0; j<width; j+=2) {

            float32 f = 0.0f;
            
            for(int k=-2; k<=+2; k++) {
                for(int l=-2; l<=+2; l++) {
                    f += X[i+k][j+l] * K[k][l];
                }
            }
            
            Y[i/2][j/2] = f;
        }
    }
}
/* ------------------------------------------------------------ */
void resample_F32(float32 **X, int height, int width, float32 **Y)
/* ------------------------------------------------------------ */
{
    int i, j;
    int ii, jj;

    float32 a,b,c,d;
    float32 aa, bb, cc, dd;

    // a . b     a u1 b
    // . . . => u2 u3 .
    // c . d     c .  .

    for(i = 0; i < height - 1; i++) {

        ii = 2 * i;

        for(j = 0; j < width - 1; j++) {

            jj = 2*j;

            a = X[i  ][j]; b = X[i  ][j+1];
            c = X[i+1][j]; d = X[i+1][j+1];

            aa = a + a;
            bb = a + b; // 2*((a+b)/2)
            cc = a + c; // 2*((a+c)/2)
            dd = (float32) (0.5f * (a + b + c+ d));

            Y[ii  ][jj] = aa; Y[ii  ][jj+1] = bb;
            Y[ii+1][jj] = cc; Y[ii+1][jj+1] = dd;
        }
        // last column

        j = width - 1; jj = 2 * j;
        
        a = X[i  ][j];
        c = X[i+1][j];
        aa = a + a;
        bb = a + c;

        Y[ii  ][jj] = aa; Y[ii  ][jj+1] = aa;
        Y[ii+1][jj] = bb; Y[ii+1][jj+1] = bb;
    }
    // last line
    i = height - 1; ii = 2 * i;
    for(j = 0; j < width - 1; j++) {
        
        jj = 2 * j;
        
        a = X[i][j  ];
        b = X[i][j+1];
        
        aa = a + a;
        bb = a + b;

        Y[ii  ][jj] = aa; Y[ii  ][jj+1] = bb;
        Y[ii+1][jj] = aa; Y[ii+1][jj+1] = bb;
    }
    // last point
    j = width - 1; jj = 2 * j;
    a = X[i][j];
    aa = a + a;

    Y[ii  ][jj] = aa; Y[ii  ][jj+1] = aa;
    Y[ii+1][jj] = aa; Y[ii+1][jj+1] = aa;
}
/* --------------------------------------------------------------------- */
void findDiff_F32(float32 **X, int h, int w, int b, float32 **Y, char *msg)
/* --------------------------------------------------------------------- */
{
    int i, j;
    float32 e = 1e-8;
    float32 x, y, d;
    for(i=0-b; i<=h-1+b; i++) {
        for(j=0-b; j<=w-1+b; j++) {
            x = X[i][j];
            y = Y[i][j];
            if(fabs(x) > 1) {
                d = fabs((x-y)/x);
            } else {
                d = fabs(x-y);
            }
            if(d > e) {
                printf("%s diff en i=%3d j=%3d : %15.10f %15.10f\n", msg, i, j, X[i][j], Y[i][j]);
            }
        }
    }
}
// -------------------------------------------------------------------------------------------------------
void convert_ui8matrix_f32matrix_norm(uint8 **X, int i0, int i1, int j0, int j1, float32 norm, float32 **Y)
// -------------------------------------------------------------------------------------------------------
{
    float32 x, y;
    float32 inv = 1.0f / norm;
    
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++) {
            x = X[i][j];
            y = inv * x;
            Y[i][j] = y;
        }
    }
}
// --------------------------------------------------------------------------------------------------------
void convert_f32matrix_ui8matrix_norm(float32 **X, int i0, int i1, int j0, int j1, float32 norm, uint8 **Y)
// --------------------------------------------------------------------------------------------------------
{
    uint8 b;
    float32 x, y;
    float32 inv = 1.0f / norm;
    
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++) {
            x = X[i][j];
            y = inv * x;         
            if(y<0)
                b = 0;
            else
                if(y>255)
                    b = 255;
            else
                b = (uint8) y;
            
            Y[i][j] = b;
        }
    }
}
// ----------------------------------------------------------------------------
void calc_L1norm(float32 **U, float32 **V, int height, int width, float32 **UV)
// ----------------------------------------------------------------------------
{
    for(int i=0; i<height; i++) {
        for(int j=0; j<width; j++) {
        
            UV[i][j] = sqrtf(U[i][j] * U[i][j] + V[i][j] * V[i][j]);
        }
    }
}