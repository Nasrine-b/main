/* ----------------------------- */
/* --- HornShunck_operator.c --- */
/* ----------------------------- */

/*
 * Copyright (c) 2017-2018, Lionel Lacassagne, all rights reserved
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
//#include "parser.h"
//#include "sequence.h"

//#include "of_config.h"
//#include "of.h"
#include "of_macro.h"
#include "of_function.h"

//#include "hs_algo.h"
#include "hs_operator.h"

#include "interp.h"
#include "interp_extra.h"

// -----------------------------------------------------------------------------------------------------------------------------------------
void hs_basic_F32(float32 **E0, float32 **E1, float32 **U0, float32 **V0, int height, int width, float32 alpha2, float32 **U1, float32 **V1)
// -----------------------------------------------------------------------------------------------------------------------------------------
{
    float32 Ex, Ey, Et; // first derivatives
    float32 u, v;       // speed (u,v) in (i,j)
    float32 um, vm;     // averages
    
    //alpha2 = 0.635;
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
            
            // horinzontal derivative
            Ex  = (E0[i][j+1] - E0[i][j]) + (E0[i+1][j+1] - E0[i+1][j]); // k
            Ex += (E1[i][j+1] - E1[i][j]) + (E1[i+1][j+1] - E1[i+1][j]); // k+1
            Ex /=  4.0f;
            
            // vertical derivative
            Ey  = (E0[i+1][j] - E0[i][j]) + (E0[i+1][j+1] - E0[i][j+1]);
            Ey += (E1[i+1][j] - E1[i][j]) + (E1[i+1][j+1] - E1[i][j+1]);
            Ey /= 4.0f;
            
            // temporal derivative
            Et  = E1[i  ][j  ] - E0[i  ][j  ];
            Et += E1[i  ][j+1] - E0[i  ][j+1];
            Et += E1[i+1][j  ] - E0[i+1][j  ];
            Et += E1[i+1][j+1] - E0[i+1][j+1];
            Et /= 4.0f;
            
            um  = (U0[i-1][j-1] + U0[i-1][j+1] + U0[i+1][j-1] + U0[i+1][j+1]) / 12.0f;
            um += (U0[i-1][j  ] + U0[i  ][j-1] + U0[i  ][j+1] + U0[i+1][j  ]) /  6.0f;
            
            vm  = (V0[i-1][j-1] + V0[i-1][j+1] + V0[i+1][j-1] + V0[i+1][j+1]) / 12.0f;
            vm += (V0[i-1][j  ] + V0[i  ][j-1] + V0[i  ][j+1] + V0[i+1][j  ]) /  6.0f;
            
            u = um - Ex * (Ex * um + Ey * vm + Et) / (alpha2 + Ex * Ex + Ey * Ey);
            v = vm - Ey * (Ex * um + Ey * vm + Et) / (alpha2 + Ex * Ex + Ey * Ey);
            
            U1[i][j] = u;
            V1[i][j] = v;
        }
    }
}
// ----------------------------------------------------------------------------------------
void hs_calc_gradientX_F32(float32 **E0, float32 **E1, int height, int width, float32 **Ex)
// ----------------------------------------------------------------------------------------
{
    float32 ex;
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
            
            // horinzontal derivative
            ex  = (E0[i][j+1] - E0[i][j]) + (E0[i+1][j+1] - E0[i+1][j]); // k
            ex += (E1[i][j+1] - E1[i][j]) + (E1[i+1][j+1] - E1[i+1][j]); // k+1
            ex /=  4.0f;
            Ex[i][j] = ex;
        }
    }
}
// ----------------------------------------------------------------------------------------
void hs_calc_gradientY_F32(float32 **E0, float32 **E1, int height, int width, float32 **Ey)
// ----------------------------------------------------------------------------------------
{
    float32 ey;
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
            
            // vertical derivative
            ey  = (E0[i+1][j] - E0[i][j]) + (E0[i+1][j+1] - E0[i][j+1]);
            ey += (E1[i+1][j] - E1[i][j]) + (E1[i+1][j+1] - E1[i][j+1]);
            ey /= 4.0f;
            Ey[i][j] = ey;
        }
    }
}
// ----------------------------------------------------------------------------------------
void hs_calc_gradientT_F32(float32 **E0, float32 **E1, int height, int width, float32 **Et)
// ----------------------------------------------------------------------------------------
{
    float32 et;
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
        
            // temporal derivative
            et  = E1[i  ][j  ] - E0[i  ][j  ];
            et += E1[i  ][j+1] - E0[i  ][j+1];
            et += E1[i+1][j  ] - E0[i+1][j  ];
            et += E1[i+1][j+1] - E0[i+1][j+1];
            et /= 4.0f;
            Et[i][j] = et;
        }
    }
}
// ----------------------------------------------------------------------------------------------------------------------
void hs_calc_gradientXYT_F32(float32 **E0, float32 **E1, int height, int width, float32 **Ex, float32 **Ey, float32 **Et)
// ----------------------------------------------------------------------------------------------------------------------
{
    float32 ex, ey, et;
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
            
            // horinzontal derivative
            ex  = (E0[i][j+1] - E0[i][j]) + (E0[i+1][j+1] - E0[i+1][j]); // k
            ex += (E1[i][j+1] - E1[i][j]) + (E1[i+1][j+1] - E1[i+1][j]); // k+1
            ex /=  4.0f;
            Ex[i][j] = ex;
            
            // vertical derivative
            ey  = (E0[i+1][j] - E0[i][j]) + (E0[i+1][j+1] - E0[i][j+1]);
            ey += (E1[i+1][j] - E1[i][j]) + (E1[i+1][j+1] - E1[i][j+1]);
            ey /= 4.0f;
            Ey[i][j] = ey;
            
            // temporal derivative
            et  = E1[i  ][j  ] - E0[i  ][j  ];
            et += E1[i  ][j+1] - E0[i  ][j+1];
            et += E1[i+1][j  ] - E0[i+1][j  ];
            et += E1[i+1][j+1] - E0[i+1][j+1];
            et /= 4.0f;
            Et[i][j] = et;
        }
    }
}
// --------------------------------------------------------------------------------------------------------------------------------------
void hs_calc_UmVm612_F32(float32** restrict U, float32** restrict V, int height, int width, float32** restrict Um, float32** restrict Vm)
// --------------------------------------------------------------------------------------------------------------------------------------
{
#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width) shared(U, V, Um, Vm) default(none)
#endif
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
        
            float32 um, vm;
            
            um  = (U[i-1][j-1] + U[i-1][j+1] + U[i+1][j-1] + U[i+1][j+1]) / 12.0f;
            um += (U[i-1][j  ] + U[i  ][j-1] + U[i  ][j+1] + U[i+1][j  ]) /  6.0f;
            
            vm  = (V[i-1][j-1] + V[i-1][j+1] + V[i+1][j-1] + V[i+1][j+1]) / 12.0f;
            vm += (V[i-1][j  ] + V[i  ][j-1] + V[i  ][j+1] + V[i+1][j  ]) /  6.0f;
            
            Um[i][j] = um;
            Vm[i][j] = vm;
        }
    }
}
// --------------------------------------------------------------------------------------------------------------------------------------
void hs_calc_UmVm12_F32(float32** restrict U, float32** restrict V, int height, int width, float32** restrict Um, float32** restrict Vm)
// --------------------------------------------------------------------------------------------------------------------------------------
{
    #ifdef OPENMP
    #pragma omp parallel for firstprivate(height, width) shared(U, V, Um, Vm) default(none)
    #endif
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
            
            float32 um, vm;
            
            um  =     U[i-1][j-1] + 2 * U[i-1][j] +        U[i-1][j+1];
            um += 2 * U[i  ][j-1] +                    2 * U[i  ][j+1];
            um +=     U[i+1][j-1] + 2 * U[i+1][j] +        U[i+1][j+1];
            um /= 12;
            
            vm  =     V[i-1][j-1] + 2 * V[i-1][j] +        V[i-1][j+1];
            vm += 2 * V[i  ][j-1] +                    2 * V[i  ][j+1];
            vm +=     V[i+1][j-1] + 2 * V[i+1][j] +        V[i+1][j+1];
            vm /= 12;
            
            Um[i][j] = um;
            Vm[i][j] = vm;
        }
    }
}
// -------------------------------------------------------------------------------------------------
void hs_calc_UmVm16_F32(float32 **U, float32 **V, int height, int width, float32 **Um, float32 **Vm)
// -------------------------------------------------------------------------------------------------
{

#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width) shared(U, V, Um, Vm) default(none)
#endif
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
        
            float32 um, vm;
            
            um  =     U[i-1][j-1] + 2 * U[i-1][j] +     U[i-1][j+1];
            um += 2 * U[i  ][j-1] + 4 * U[i  ][j] + 2 * U[i  ][j+1];
            um +=     U[i+1][j-1] + 2 * U[i+1][j] +     U[i+1][j+1];
            
            vm  =     V[i-1][j-1] + 2 * V[i-1][j] +     V[i-1][j+1];
            vm += 2 * V[i  ][j-1] + 4 * V[i  ][j] + 2 * V[i  ][j+1];
            vm +=     V[i+1][j-1] + 2 * V[i+1][j] +     V[i+1][j+1];

            um /= 16.0f;
            vm /= 16.0f;

            Um[i][j] = um;
            Vm[i][j] = vm;
        }
    }
}
// -----------------------------------------------------------------------------------------------
void hs_calc_UmVm_F32(float32 **U, float32 **V, int height, int width, float32 **Um, float32 **Vm)
// -----------------------------------------------------------------------------------------------
{
    //hs_calc_UmVm16_F32(U, V, height, iwdth, Um, Vm); // fast
    hs_calc_UmVm12_F32(U, V, height, width, Um, Vm); // accurate
}

// ----------------------------------------------------------------------------------------------------------------------------
void hs_calc_UV0_F32(float32 **Ex, float32 **Ey, float32 **Et, int height, int width, float32 alpha2, float32 **U, float32 **V)
// ----------------------------------------------------------------------------------------------------------------------------
{
    // special case for t=0
    
#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width, alpha2) shared(Ex, Ey, Et, U, V) default(none)
#endif

    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
            
            float32 ex  = Ex[i][j];
            float32 ey  = Ey[i][j];
            float32 et  = Et[i][j];
            
            // um = 0, vm = 0
            float32 u = - (ex * et) / (alpha2 + ex * ex + ey * ey);
            float32 v = - (ey * et) / (alpha2 + ex * ex + ey * ey);
            
            U[i][j] = u;
            V[i][j] = v;
        }
    }
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------
void hs_calc_UV_F32(float32 **Ex, float32 **Ey, float32 **Et, float32 **Um, float32 **Vm, int height, int width, float32 alpha2, float32 **U, float32 **V)
// -------------------------------------------------------------------------------------------------------------------------------------------------------
{
    // general case for t
   
#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width, alpha2) shared(Ex, Ey, Et, Um, Vm, U, V) default(none)
#endif
   
    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
            
            float32 ex  = Ex[i][j];
            float32 ey  = Ey[i][j];
            float32 et  = Et[i][j];
            
            float32 um  = Um[i][j];
            float32 vm  = Vm[i][j];
            
            float32 u = um - ex * (ex * um + ey * vm + et) / (alpha2 + ex * ex + ey * ey);
            float32 v = vm - ey * (ex * um + ey * vm + et) / (alpha2 + ex * ex + ey * ey);
            
            U[i][j] = u;
            V[i][j] = v;
        }
    }
}
// ----------------------------------------------------------------------------------------------------------------------------------
void hs_calc_UV_direct_F32(float32 **Ex, float32 **Ey, float32 **Et, int height, int width, float32 alpha2, float32 **U, float32 **V)
// ----------------------------------------------------------------------------------------------------------------------------------
{
    // sans passer par Um,Vm

#ifdef OPENMP
#pragma omp parallel for firstprivate(height, width, alpha2) shared(Ex, Ey, Et, U, V) default(none)
#endif

    for(int i = 0; i < height; i++) {
        for(int j = 0; j< width; j++) {
            
            float32 ex  = Ex[i][j];
            float32 ey  = Ey[i][j];
            float32 et  = Et[i][j];
            
            /*
             um  =        U[i-1][j-1] + 2.0f * U[i-1][j] +        U[i-1][j+1];
             um += 2.0f * U[i  ][j-1] + 4.0f * U[i  ][j] + 2.0f * U[i  ][j+1];
             um +=        U[i+1][j-1] + 2.0f * U[i+1][j] +        U[i+1][j+1];
             um /= 16.0f;
             
             vm  =        V[i-1][j-1] + 2.0f * V[i-1][j] +        V[i-1][j+1];
             vm += 2.0f * V[i  ][j-1] + 4.0f * V[i  ][j] + 2.0f * V[i  ][j+1];
             vm +=        V[i+1][j-1] + 2.0f * V[i+1][j] +        V[i+1][j+1];
             vm /= 16.0f;*/
            
            
            float32 um  = U[i-1][j-1] + U[i-1][j+1] + U[i+1][j-1] + U[i+1][j+1];
            um += 2.0f * (U[i-1][j  ] + U[i  ][j-1] + U[i  ][j+1] + U[i+1][j  ]);
            um /= 12.0f;
            
            float32 vm  = V[i-1][j-1] + V[i-1][j+1] + V[i+1][j-1] + V[i+1][j+1];
            vm += 2.0f * (V[i-1][j  ] + V[i  ][j-1] + V[i  ][j+1] + V[i+1][j  ]);
            vm /= 12.0f;/**/
            
            float32 u = um - ex * (ex * um + ey * vm + et) / (alpha2 + ex * ex + ey * ey);
            float32 v = vm - ey * (ex * um + ey * vm + et) / (alpha2 + ex * ex + ey * ey);
            
            U[i][j] = u;
            V[i][j] = v;
        }
    }
}
// ---------------------------------------------------------------------------------------------------
void hs_accumulate_UV_F32(float32 **dU, float32 **dV, int height, int width, float32 **U, float32 **V)
// ---------------------------------------------------------------------------------------------------
{
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
        
            U[i][j] += dU[i][j];
        }
    }
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            
            V[i][j] += dV[i][j];
        }
    }
}
// ------------------------------------------------------------------------------------------------------------------------
void hs_add_UV_F32(float32 **U, float32 **V, float32 **dU, float32 **dV, int height, int width, float32 **Ut, float32 **Vt)
// ------------------------------------------------------------------------------------------------------------------------
{
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            
            Ut[i][j] = U[i][j] + dU[i][j];
        }
    }
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
        
            Vt[i][j] = V[i][j] + dV[i][j];
        }
    }
}
// -----------------------------------------------------------------
void hs_clip_UV_F32(float32 **U, float32 **V, int height, int width)
// -----------------------------------------------------------------
{
    float32 di0, di1;
    float32 dj0, dj1;
    
    //return;
    
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            
            dj0 = (float32) (0.0f  - j);
            dj1 = (float32) (width - j);
            
            if(U[i][j] < dj0) U[i][j] = dj0;
            if(U[i][j] > dj1) U[i][j] = dj1;
            
            /*if(U[i][j] < dj0) {
                printf("U[%2d][%2d] = %8.3f <- %8.3f\n", i, j, U[i][j], dj0);
                U[i][j] = dj0;
            }
            if(U[i][j] > dj1) {
                printf("U[%2d][%2d] = %8.3f <- %8.3f\n", i, j, U[i][j], dj1);
                U[i][j] = dj1;
            }*/
            
            //printf("%8.3f %8.3f ", dj0, dj1);
        }
        //puts("");
    }
    //puts("");
    //exit(-1);
    
    //puts("V processing");
    for(int i = 0; i < height; i++) {
        
        di0 = (float32) (0.0f   - i);
        di1 = (float32) (height - i);
        
        for(int j = 0; j < width; j++) {

            if(V[i][j] < di0) V[i][j] = di0;
            if(V[i][j] > di1) V[i][j] = di1;
            
            /*if(V[i][j] < di0) {
                printf("V[%2d][%2d] = %8.3f <- %8.3f\n", i, j, V[i][j], di0);
                V[i][j] = di0;
            }
            if(V[i][j] > di1) {
                printf("V[%2d][%2d] = %8.3f <- %8.3f\n", i, j, V[i][j], di1);
                V[i][j] = di1;
            }*/
        }
    }
    //exit(-1);
}
// -------------------------------------------------------------------------------------------------------------------------------------------------------
void hs_step0_F32(float32 **E0, float32 **E1, int height, int width, float32 alpha2, float32 **Ex, float32 **Ey, float32 **Et, float32 **U1, float32 **V1)
// -------------------------------------------------------------------------------------------------------------------------------------------------------
{
    // only one time, before the iterations
    
    //hs_calc_gradientX_F32(E0, E1, height, width, Ex);
    //hs_calc_gradientY_F32(E0, E1, height, width, Ey);
    //hs_calc_gradientT_F32(E0, E1, height, width, Et);
    hs_calc_gradientXYT_F32(E0, E1, height, width, Ex, Ey, Et);
    
    // no need to compute UmVm because U(t=0)=0 and V(t=0)=0
    hs_calc_UV0_F32(Ex, Ey, Et, height, width, alpha2, U1, V1);
}
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void hs_1step_F32(float32 **U0, float32 **V0, float32 **Ex, float32 **Ey, float32 **Et, int height, int width, float32 alpha2, float32 **Um, float32 **Vm, float32 **U1, float32 **V1)
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
{
    hs_calc_UmVm_F32(U0, V0, height, width, Um, Vm);
    hs_calc_UV_F32(Ex, Ey, Et, Um, Vm, height, width, alpha2, U1, V1);
}
// -------------------------
//void HS_interpLinear(OF *of)
// -------------------------
/*{
    int level = of->level;
    int h   = of->height >> level;
    int w   = of->width  >> level;
    
    float32 **U = of->U[level];
    float32 **V = of->V[level];
    
    uint8 **I, **Irec;
    float32 **F, **Frec;
    
    if(of->iter <= of->nb_iter) {
        // les vecteurs temporaires (du+u,dv+v) sont dans (ua,va) (pour debug)
        U = of->Ut[level];
        V = of->Vt[level];
    } else {
        // les vecteurs vitesses definitifs (du+u,dv+v) sont accumules dans (u,v)
        U = of->U[level];
        V = of->V[level];
    }
    
    if(of->calc_type == CALC_I8) {
        I    = of->I2[level];
        Irec = of->I2rec[level];
        printf("[HS_interpLinear]: interpLinear_I8 n'existe plus\n"); exit(-1);
        //interpLinear_I8(I, h, w, U, V, Irec);
    }
    
    if(of->calc_type == CALC_F32) {
        F    = of->F2[level];
        Frec = of->F2rec[level];
        OF_EXTEND(extendDataD2_F32(F, h, w));
        interpLinear_F32(F, U, V, h, w, Frec);
        
    }
}*/
// ------------------------
//void HS_interpCubic(OF *of)
// ------------------------
/*{
    int level    = of->level;
    int iter     = of->iter;
    int nb_iter  = of->nb_iter;
    int h        = of->height >> level;
    int w        = of->width  >> level;
    
    float32 **U = of->U[level];
    float32 **V = of->V[level];
    
    uint8   **I, **Irec;
    float32 **F, **Frec;
    
    if( (iter >= 1) && (iter <= nb_iter)) {
        // cas principal = dans la boucle
        // les vecteurs temporaires (du+u,dv+v) sont dans (ua,va) (pour debug)
        U = of->Ut[level];
        V = of->Vt[level];
    } else {
        // cas particuliers prologue ou epilogue
        // les vecteurs vitesses definitifs (du+u,dv+v) sont accumules dans (u,v)
        U = of->U[level];
        V = of->V[level];
    }
    
    if(of->calc_type == CALC_I8) {
        I    = of->I2[level];
        Irec = of->I2rec[level];
        printf("[HS_interpCubic]: interpCubic_I8 n'existe plus"); exit(-1);
        //interpCubic_I8(I, h, w, U, V, Irec);
    }
    
    if(of->calc_type == CALC_F32) {
        F    = of->F2[level];
        Frec = of->F2rec[level];
        OF_EXTEND(extendDataD2_F32(F, h, w));
        
#ifdef FOLKI_ENABLE_SIMD_INTERP
        interpCubicLUT_vF32(of->vF2[level], h, w, U, V, of->interpCubicLUT_vF32, of->interpSizeLUT, of->F2rec[level]);
#else
        interpCubic_F32(F, U, V, h, w, Frec);
#endif
    }
}*/
