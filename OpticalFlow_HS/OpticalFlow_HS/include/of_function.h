/* --------------------- */
/* --- OF_function.h --- */
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

#ifndef __OF_FUNCTION_H__
#define __OF_FUNCTION_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

int is_nan_F32(float32 f);
int is_nan_F64(float64 f);
    int isnan_f32matrix(float32 **m, int i0, int i1, int j0, int j1, char *name);
    
void makeBorder_I8(uint8 **T, int i0, int i1, int j0, int j1, int border);
void zeroBorder_I8(uint8 **T, int i0, int i1, int j0, int j1, int border);
    
void extendDataZ2_I8(uint8 **X, int h, int w);
void extendDataD2_I8(uint8 **X, int h, int w);
    
void makeBorder_F32(float32 **T, int i0, int i1, int j0, int j1, int border);
void zeroBorder_F32(float32 **T, int i0, int i1, int j0, int j1, int border);

void extendDataZ2_F32(float32 **X, int h, int w);
void extendDataD2_F32(float32 **X, int h, int w);

float32 max_matrix_pos(float32 **X, int i0, int i1, int j0, int j1, int *ipos, int *jpos);

void display_ui16matrix_quantif(uint16 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name);
void display_si16matrix_quantif(sint16 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name);
void display_si32matrix_quantif(sint32 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name);

void convert_f32matrix_ui8matrix_sat(float32 **X, int i0, int i1, int j0, int j1, uint8 **Y);

// --------------------------- //
// -- kernel initialisation -- //
// --------------------------- //

void initKernel_F32(float32 **K, float32 c0, float32 c1, float32 c2, float32 c3, float32 c4);
void initBurtKernel_F32    (float32 a, float32 **K);
void initGauss5Kernel_F32  (float32 **K);
void initGauss3Kernel_F32  (float32 **K);
void initAverage5Kernel_F32(float32 **K);
void initAverage3Kernel_F32(float32 **K);

// --------------------------- //
// --- fonctions de calcul --- //
// --------------------------- //

void convolution5x5_F32(float32 **X, int height, int width, float32 **K,  float32 **Y);
void convolution5x5_decimation_I8_F32_I8   (uint8   **X, int height, int width, float32 **K,  uint8    **Y);
void convolution5x5_decimation_F32_F32_F32(float32 **X, int height, int width, float32 **K,  float32 **Y);

void resample_F32(float32 **X, int height, int width, float32 **Y);

float32 calcMSE_F32(float32 **X, int height, int width, float32 **Y);
float32 calcMSE_I8 (uint8    **X, int height, int width, uint8    **Y);
double calcAngle(double x0, double y0, double x1, double y1);
float32 calcAngle_F32(float32 **X0, float32 **Y, int height, int width, float32 **X1, float32 **Y1);

// ------------------ //
// --- conversion --- //
// ------------------ //
void convert_ui8matrix_f32matrix_norm(uint8   **X, int i0, int i1, int j0, int j1, float32 norm, float32 **Y);
void convert_f32matrix_ui8matrix_norm(float32 **X, int i0, int i1, int j0, int j1, float32 norm, uint8   **Y);

void convert_f32matrix_ui8matrix_sat(float32 **X, int i0, int i1, int j0, int j1, uint8 **Y);

// ------------ //
// --- test --- //
// ------------ //
void findDiff_F32(float32 **X, int h, int w, int b, float32 **Y, char *msg);

void calc_L1norm(float32 **U, float32 **V, int height, int width, float32 **UV);

#ifdef __cplusplus
}
#endif

#endif
