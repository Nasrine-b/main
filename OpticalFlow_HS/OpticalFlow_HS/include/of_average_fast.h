/* ------------------------- */
/* --- of_average_fast.h --- */
/* ------------------------- */

/*
 * Lukas Kanade main header
 */

/*
 * Copyright (c) 2006-2008 Lionel Lacassagne, IEF, all rights reserved
 */

/*
 * code CONFIDENTIEL LABORATOIRE
 * Export strictement interdit en dehors de l'IEF
 * tout contrevenant sera poursuivi
 */

#ifndef __OF_AVERAGE_FAST_H__
#define __OF_AVERAGE_FAST_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

#define sumH1(X, i, j) X[i][j-1] + X[i][j] + X[i][j+1]
#define sumH2(X, i, j) X[i][j-2] + sumH1(X,i,j) + X[i][j+2]
#define sumH3(X, i, j) X[i][j-3] + sumH2(X,i,j) + X[i][j+3]
#define sumH4(X, i, j) X[i][j-4] + sumH3(X,i,j) + X[i][j+4]
#define sumH5(X, i, j) X[i][j-5] + sumH4(X,i,j) + X[i][j+5]
#define sumH6(X, i, j) X[i][j-6] + sumH5(X,i,j) + X[i][j+6]
#define sumH7(X, i, j) X[i][j-7] + sumH6(X,i,j) + X[i][j+7]
#define sumH8(X, i, j) X[i][j-8] + sumH7(X,i,j) + X[i][j+8]
#define sumH9(X, i, j) X[i][j-9] + sumH8(X,i,j) + X[i][j+9]
#define sumH10(X, i, j) X[i][j-10] + sumH9(X,i,j) + X[i][j+10]
    
#define sumV1(X,i,j) X[i-1][j] + X[i][j] + X[i+1][j]
#define sumV2(X,i,j) X[i-2][j] + sumV1(X,i,j) + X[i+2][j]
#define sumV3(X,i,j) X[i-3][j] + sumV2(X,i,j) + X[i+3][j]
#define sumV4(X,i,j) X[i-4][j] + sumV3(X,i,j) + X[i+4][j]
#define sumV5(X,i,j) X[i-5][j] + sumV4(X,i,j) + X[i+5][j]
#define sumV6(X,i,j) X[i-6][j] + sumV5(X,i,j) + X[i+6][j]
#define sumV7(X,i,j) X[i-7][j] + sumV6(X,i,j) + X[i+7][j]
#define sumV8(X,i,j) X[i-8][j] + sumV7(X,i,j) + X[i+8][j]
#define sumV9(X,i,j) X[i-9][j] + sumV8(X,i,j) + X[i+9][j]
    
#define calc_factorH1(X, i, j) X[i][j] + X[i][j+1]
#define calc_factorH2(X, i, j) X[i][j-1] + calc_factorH1(X, i, j) + X[i][j+2]
#define calc_factorH3(X, i, j) X[i][j-2] + calc_factorH2(X, i, j) + X[i][j+3]
#define calc_factorH4(X, i, j) X[i][j-3] + calc_factorH3(X, i, j) + X[i][j+4]
#define calc_factorH5(X, i, j) X[i][j-4] + calc_factorH4(X, i, j) + X[i][j+5]
#define calc_factorH6(X, i, j) X[i][j-5] + calc_factorH5(X, i, j) + X[i][j+6]
#define calc_factorH7(X, i, j) X[i][j-6] + calc_factorH6(X, i, j) + X[i][j+7]
#define calc_factorH8(X, i, j) X[i][j-7] + calc_factorH7(X, i, j) + X[i][j+8]
#define calc_factorH9(X, i, j) X[i][j-8] + calc_factorH8(X, i, j) + X[i][j+9]
    
#define sumH1_LU2(X, i, j, s0, s1) factor = calc_factorH1(X,i,j); s0 = factor + X[i][j-1]; s1 = factor + X[i][j+2]
#define sumH2_LU2(X, i, j, s0, s1) factor = calc_factorH2(X,i,j); s0 = factor + X[i][j-2]; s1 = factor + X[i][j+3]
#define sumH3_LU2(X, i, j, s0, s1) factor = calc_factorH3(X,i,j); s0 = factor + X[i][j-3]; s1 = factor + X[i][j+4]
#define sumH4_LU2(X, i, j, s0, s1) factor = calc_factorH4(X,i,j); s0 = factor + X[i][j-4]; s1 = factor + X[i][j+5]
#define sumH5_LU2(X, i, j, s0, s1) factor = calc_factorH5(X,i,j); s0 = factor + X[i][j-5]; s1 = factor + X[i][j+6]
#define sumH6_LU2(X, i, j, s0, s1) factor = calc_factorH6(X,i,j); s0 = factor + X[i][j-6]; s1 = factor + X[i][j+7]
#define sumH7_LU2(X, i, j, s0, s1) factor = calc_factorH7(X,i,j); s0 = factor + X[i][j-7]; s1 = factor + X[i][j+8]
#define sumH8_LU2(X, i, j, s0, s1) factor = calc_factorH8(X,i,j); s0 = factor + X[i][j-8]; s1 = factor + X[i][j+9]
#define sumH9_LU2(X, i, j, s0, s1) factor = calc_factorH9(X,i,j); s0 = factor + X[i][j-9]; s1 = factor + X[i][j+10]
    
#define calc_factorV1(X, i, j) X[i][j] + X[i+1][j]
#define calc_factorV2(X, i, j) X[i-1][j] + calc_factorV1(X, i, j) + X[i+2][j]
#define calc_factorV3(X, i, j) X[i-2][j] + calc_factorV2(X, i, j) + X[i+3][j]
#define calc_factorV4(X, i, j) X[i-3][j] + calc_factorV3(X, i, j) + X[i+4][j]
#define calc_factorV5(X, i, j) X[i-4][j] + calc_factorV4(X, i, j) + X[i+5][j]
#define calc_factorV6(X, i, j) X[i-5][j] + calc_factorV5(X, i, j) + X[i+6][j]
#define calc_factorV7(X, i, j) X[i-6][j] + calc_factorV6(X, i, j) + X[i+7][j]
#define calc_factorV8(X, i, j) X[i-7][j] + calc_factorV7(X, i, j) + X[i+8][j]
#define calc_factorV9(X, i, j) X[i-8][j] + calc_factorV8(X, i, j) + X[i+9][j]
    
#define sumV1_LU2(X, i, j, s0, s1) factor = calc_factorV1(X,i,j); s0 = factor + X[i-1][j]; s1 = factor + X[i+2][j]
#define sumV2_LU2(X, i, j, s0, s1) factor = calc_factorV2(X,i,j); s0 = factor + X[i-2][j]; s1 = factor + X[i+3][j]
#define sumV3_LU2(X, i, j, s0, s1) factor = calc_factorV3(X,i,j); s0 = factor + X[i-3][j]; s1 = factor + X[i+4][j]
#define sumV4_LU2(X, i, j, s0, s1) factor = calc_factorV4(X,i,j); s0 = factor + X[i-4][j]; s1 = factor + X[i+5][j]
#define sumV5_LU2(X, i, j, s0, s1) factor = calc_factorV5(X,i,j); s0 = factor + X[i-5][j]; s1 = factor + X[i+6][j]
#define sumV6_LU2(X, i, j, s0, s1) factor = calc_factorV6(X,i,j); s0 = factor + X[i-6][j]; s1 = factor + X[i+7][j]
#define sumV7_LU2(X, i, j, s0, s1) factor = calc_factorV7(X,i,j); s0 = factor + X[i-7][j]; s1 = factor + X[i+8][j]
#define sumV8_LU2(X, i, j, s0, s1) factor = calc_factorV8(X,i,j); s0 = factor + X[i-8][j]; s1 = factor + X[i+9][j]
#define sumV9_LU2(X, i, j, s0, s1) factor = calc_factorV9(X,i,j); s0 = factor + X[i-9][j]; s1 = factor + X[i+10][j]
    
    /*
     #define hSum4_2D(i,j) X[i][j] + X[i][j+1] + X[i][j+2] + X[i][j+3]
     #define hSum8_2D(i,j) X[i][j] + X[i][j+1] + X[i][j+2] + X[i][j+3] + X[i][j+4] + X[i][j+5] + X[i][j+6] + X[i][j+7]
     #define hSum10_2D(i,j) X[i][j] + X[i][j+1] + X[i][j+2] + X[i][j+3] + X[i][j+4] + X[i][j+5] + X[i][j+6] + X[i][j+7] + X[i][j+8] + X[i][j+9]
     
     #define vhSum4_2D(i,j) hSum4_2D(i,j) + hSum4_2D(i+1,j) + hSum4_2D(i+2,j) + hSum4_2D(i+3,j)
     #define vhSum8_2D(i,j) hSum4_2D(i,j) + hSum4_2D(i+1,j) + hSum4_2D(i+2,j) + hSum4_2D(i+3,j) + hSum4_2D(i+4,j) + hSum4_2D(i+5,j) + hSum4_2D(i+6,j) + hSum4_2D(i+7,j)
     #define vhSum10_2D(i,j) hSum4_2D(i,j) + hSum4_2D(i+1,j) + hSum4_2D(i+2,j) + hSum4_2D(i+3,j) + hSum4_2D(i+4,j) + hSum4_2D(i+5,j) + hSum4_2D(i+6,j) + hSum4_2D(i+7,j)+ hSum4_2D(i+8,j) + hSum4_2D(i+9,j)
     
     #define hSum4_1D(i,j) X[i][j] + X[i][j+1] + X[i][j+2] + X[i][j+3]
     #define hSum8_1D(i,j) X[i][j] + X[i][j+1] + X[i][j+2] + X[i][j+3] + X[i][j+4] + X[i][j+5] + X[i][j+6] + X[i][j+7]
     #define hSum10_1D(i,j) X[i][j] + X[i][j+1] + X[i][j+2] + X[i][j+3] + X[i][j+4] + X[i][j+5] + X[i][j+6] + X[i][j+7] + X[i][j+8] + X[i][j+9]
     
     #define vhSum4_1D(i,j)  hSum4_1D(i,j) + hSum4_1D(i+1,j) + hSum4_1D(i+2,j) + hSum4_1D(i+3,j)
     #define vhSum8_1D(i,j)  hSum4_1D(i,j) + hSum4_1D(i+1,j) + hSum4_1D(i+2,j) + hSum4_1D(i+3,j) + hSum4_1D(i+4,j) + hSum4_1D(i+5,j) + hSum4_1D(i+6,j) + hSum4_1D(i+7,j)
     #define vhSum10_1D(i,j) hSum4_1D(i,j) + hSum4_1D(i+1,j) + hSum4_1D(i+2,j) + hSum4_1D(i+3,j) + hSum4_1D(i+4,j) + hSum4_1D(i+5,j) + hSum4_1D(i+6,j) + hSum4_1D(i+7,j)+ hSum4_1D(i+8,j) + hSum4_1D(i+9,j)
     */
    
//void moyenneur_Basic_F32(float32 **X, int height, int width, int rayon, float32 **Y);
//void moyenneur_Opt_F32(float32 **X, int height, int width, int rayon, float32 **Y);
    
// fonction principale
void moyenneur_F32  (float32 **X, int height, int width, int rayon, float32 **S, float32 **Y);   
void somme_Basic_F32(float32 **X, int height, int width, int rayon, float32 **Y);

// -------------------
// fenetres glissantes
// -------------------

// 1 passe
void slidingAverage_F32           (float32 **X, int h, int w, int r, float32 **S);

// 2 passe
void slidingSum_F32               (float32 **X, int h, int w, int r, float32 **S);
void averageSlidingSum_F32        (float32 **S, int h, int w, int r, float32 **A); // normalisation
void centeredAverageSlidingSum_F32(float32 **X, int h, int w, int r, float32 **S, float32 **A); // fonction principale

void slidingSum2_F32               (float32 **X, int h, int w, int r, float32 **S1, float32 **S2);
void centeredAverageSlidingSum2_F32(float32 **X, int h, int w, int r, float32 **S, float32 **Y); // fonction principale

void slidingSum_F64               (float32 **X, int h, int w, int r, float64 **S);
void averageSlidingSum_F64        (float64 **S, int h, int w, int r, float32 **A);
void centeredAverageSlidingSum_F64(float32 **X, int h, int w, int r, float64 **S, float32 **A);  // fonction principale

// -------------------------------
// moyennage par filtre 1D et macro
// --------------------------------

void average1M_F32    (float32 **X, int h, int w, int b, float32 **S, float32 **Y); // macro
void average1MU2_F32  (float32 **X, int h, int w, int b, float32 **S, float32 **Y); // macro + Loop Unroll 2
void average1MU2X2_F32(float32 **X, int h, int w, int b, float32 **S, float32 **Y); // macro + Loop Unroll 2x2
void average1MP_F32   (float32 **X, int h, int w, int b, float32 **S, float32 **Y); // macro + Pipeline
    
// ---------------
// Sommes cumulees (images integrales)
// ---------------

// allocations sur [0-(b+1).. h-1+b][0-(b+1).. w-1+b]
void cumulatedSum_F32               (float32** restrict X, int h, int w, int b, float32** restrict S); // image integrale
void averageCumulatedSum_F32        (float32** restrict S, int h, int w, int b, float32** restrict A); // normalisation
void centeredAverageCumulatedSum_F32(float32** restrict X, int h, int w, int b, float32** restrict S, float32** restrict A); // fonction principale

void cumulatedSum_F64               (float32** restrict X, int h, int w, int b, float64** restrict S); // image integrale
void averageCumulatedSum_F64        (float64** restrict S, int h, int w, int b, float32** restrict A); // normalisation
void centeredAverageCumulatedSum_F64(float32** restrict X, int h, int w, int b, float64** restrict S, float32** restrict A);  // fonction principale

#ifdef __cplusplus
}
#endif

#endif
