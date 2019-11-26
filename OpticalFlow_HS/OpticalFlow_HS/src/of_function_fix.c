/* ------------------------- */
/* --- OF_function_fix.c --- */
/* ------------------------- */

/*
 * fonctions generalistes d'OF en virgule fixe
 */

/*
 * Copyright (c) 2018-2018 Lionel Lacassagne, all rights reserved, LIP6, SU, CNRS
 */

/*
 * 2018-12-11 creation
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
#include "of_function_fix.h"

// -------------------------------------------------------------------------------------------------------------
void display_ui16matrix_quantif(uint16 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name)
// -------------------------------------------------------------------------------------------------------------
{
    if(name != NULL) puts(name);
    
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++) {
            
            uint16 x = X[i][j];
            float32 f = (float32) x;
            printf(format, f / Q);
        }
        putchar('\n');
    }
    putchar('\n');
}
// -------------------------------------------------------------------------------------------------------------
void display_si16matrix_quantif(sint16 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name)
// -------------------------------------------------------------------------------------------------------------
{
    if(name != NULL) puts(name);
    
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++) {
            
            sint16 x = X[i][j];
            float32 f = (float32) x;
            printf(format, f / Q);
        }
        putchar('\n');
    }
    putchar('\n');
}
// -------------------------------------------------------------------------------------------------------------
void display_si32matrix_quantif(sint32 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name)
// -------------------------------------------------------------------------------------------------------------
{
    if(name != NULL) puts(name);
    
    for(int i=i0; i<=i1; i++) {
        for(int j=j0; j<=j1; j++) {
        
            sint32 x = X[i][j];
            float32 f = (float32) x;
            printf(format, f / Q);
        }
        putchar('\n');
    }
    putchar('\n');
}
// --------------------------------------------------------------------------------------------------------
void convert_ui8matrix_f32matrix_quantif(uint8 **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y)
// --------------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            Y[i][j] = X[i][j] / Q;
        }
    }
}
// ----------------------------------------------------------------------------------------------------------
void convert_ui16matrix_f32matrix_quantif(uint16 **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y)
// ----------------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            Y[i][j] = X[i][j] / Q;
        }
    }
}
// ----------------------------------------------------------------------------------------------------------
void convert_si16matrix_f32matrix_quantif(sint16 **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y)
// ----------------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            Y[i][j] = X[i][j] / Q;
        }
    }
}
// ----------------------------------------------------------------------------------------------------------
void convert_si32matrix_f32matrix_quantif(sint32 **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y)
// ----------------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            Y[i][j] = X[i][j] / Q;
        }
    }
}
// ----------------------------------------------------------------------------------------------------------
void convert_si64matrix_f32matrix_quantif(sint64 **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y)
// ----------------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            Y[i][j] = X[i][j] / Q;
        }
    }
}
// --------------------------------------------------------------------------------------------------------
void convert_f32matrix_ui8matrix_quantif(float32 **X, int i0, int i1, int j0, int j1, float32 Q, uint8 **Y)
// --------------------------------------------------------------------------------------------------------
{
    for(int i = i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            if     (X[i][j] < 0.0f) Y[i][j] = 0;
            else if(X[i][j] > 1.0f) Y[i][j] = (uint8) Q;
            else                    Y[i][j] = (uint8) lroundf(X[i][j] * Q);
        }
    }
}
// -----------------------------------------------------------------------------------------
void convert_ui16matrix_ui8matrix_sat(uint16 **X, int i0, int i1, int j0, int j1, uint8 **Y)
// -----------------------------------------------------------------------------------------
{
    sint16 x;
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
// -------------------------------------------------------
void quantifKernel5_q16(float32 **Kf, int q, uint16 **Kq)
// ------------------------------------------------------
{
    // quantification et sommation a la volee
    
    int Q = 1 << q;
    float32 f, x;
    uint16 y, s = 0;
    
    for(int i = -2; i<= +2; i++) {
        for(int j = -2; j<= +2; j++) {
            f = Kf[i][j];
            x = f * Q;
            y = roundf(x);
            s += y;
            Kq[i][j] = y;
        }
    }
    //si la somme depasse Q, l'exces est retranche au coefficient central
    Kq[0][0] -= (s - Q);
}
// ---------------------------------------------------------------------------------------------------------
void convolution5x5_decimation_U16_I64_U16(uint16 **X, int height, int width, uint16 **K, int q, uint16 **Y)
// ---------------------------------------------------------------------------------------------------------
{
    // arrondi
    sint64 r = 1 << (q-1);
    
    #ifdef OPENMP
    #pragma omp parallel for (shared X, K, Y) firstprivate(height, width, q, r) default(none)
    #endif

    for(int i=0; i<height; i+=2) {
        for(int j=0; j<width; j+=2) {

            sint64 s = 0;
            
            for(int k=-2; k<=+2; k++) {
                for(int l=-2; l<=+2; l++) {
                    s += X[i+k][j+l] * K[k][l];
                }
            }
            s += r; // arrondi simple car toujours positif
            
            Y[i/2][j/2] = (uint16) (s >> q);
        }
    }
}

// -------------------------------------------------------------
void resample_I32(sint32 **X, int height, int width, sint32 **Y)
// -------------------------------------------------------------
{
    int i, j;
    int ii, jj;

    sint32 a,b,c,d;
    sint32 aa, bb, cc, dd;

    // a . b     a u1 b
    // . . . => u2 u3 .
    // c . d     c .  .
#ifdef OPENMP
#pragma omp parallel for (shared X, Y) firstprivate(height, width, i, j, ii, jj, a, b, c, d, aa, bb, cc, dd) default(none)
#endif

    for(i = 0; i < height - 1; i++) {

        ii = 2 * i;

        for(j = 0; j < width - 1; j++) {

            jj = 2*j;

            a = X[i  ][j]; b = X[i  ][j+1];
            c = X[i+1][j]; d = X[i+1][j+1];

            aa = a + a;
            bb = a + b; // 2*((a+b)/2)
            cc = a + c; // 2*((a+c)/2)
            dd = (a + b + c+ d + 1) / 2;

            Y[ii  ][jj] = aa; Y[ii  ][jj+1] = bb;
            Y[ii+1][jj] = cc; Y[ii+1][jj+1] = dd;
        }
        
        // last column
        j = width - 1; jj = 2 * j;
        
        a = X[i  ][j];
        c = X[i+1][j];
        aa = a + a;
        bb = a + c;

        // jj+1 = 2(width-1)+1 = 2 width - 1 = OK
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
// ------------------------------------------------------------------------
void makeBorder_I16(uint16 **T, int i0, int i1, int j0, int j1, int border)
// ------------------------------------------------------------------------
{
    // --------------- //
    // -- bord haut -- //
    // --------------- //
    
    for(int k = 1; k <= border; k++) {
        for(int j = j0; j <= j1; j++) {
            T[i0-k][j] = T[i0][j];
        }
    }
    
    // -------------- //
    // -- bord bas -- //
    // -------------- //
    
    for(int k = 1; k <= border; k++) {
        for(int j = j0; j <= j1; j++) {
            T[i1+k][j] = T[i1][j];
        }
    }
    
    // --------------------------- //
    // -- bords gauche et droit -- //
    // --------------------------- //
    
    for(int i = i0-border; i <= i1+border; i++) {
        for(int k = 1; k <= border; k++) {
            T[i][j0-k] = T[i][j0];
            T[i][j1+k] = T[i][j1];
        }
    }
}

// --------------------------------------------
void extendDataD2_I16(uint16 **X, int h, int w)
// --------------------------------------------
{
    // copie de extendDataD2_F32
    
    // ERREUR DE CONCEPTION
    // R = 2
    // sauf si pour F uniquement (ie pour interpolation bi-cubique necessitant un rayon de 2)
    makeBorder_I16(X, 0, h-1, 0, w-1, 2);
}
