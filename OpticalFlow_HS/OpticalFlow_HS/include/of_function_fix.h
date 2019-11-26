/* --------------------- */
/* --- OF_function.h --- */
/* --------------------- */

/*
 * fonctions generalistes d'OF en virgule fixe
 */

/*
 * Copyright (c) 2018-2018 Lionel Lacassagne, all rights reserved, LIP6, SU, CNRS
 */

/*
 * 2018-12-11 creation
 */

#ifndef __OF_FUNCTION_FIX_H__
#define __OF_FUNCTION_FIX_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

void display_ui16matrix_quantif(uint16 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name);
void display_si16matrix_quantif(sint16 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name);
void display_si32matrix_quantif(sint32 **X, int i0, int i1, int j0, int j1, float32 Q, char *format, char *name);

void convert_ui8matrix_f32matrix_quantif (uint8   **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y);
void convert_ui16matrix_f32matrix_quantif(uint16  **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y);
void convert_si16matrix_f32matrix_quantif(sint16  **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y);
void convert_si32matrix_f32matrix_quantif(sint32  **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y);
void convert_si64matrix_f32matrix_quantif(sint64  **X, int i0, int i1, int j0, int j1, float32 Q, float32 **Y);
void convert_f32matrix_ui8matrix_quantif (float32 **X, int i0, int i1, int j0, int j1, float32 Q, uint8   **Y);

void convert_ui16matrix_ui8matrix_sat(uint16 **X, int i0, int i1, int j0, int j1, uint8 **Y);

// --------------------------- //
// -- kernel initialisation -- //
// --------------------------- //

void quantifKernel5_q16(float32 **Kf, int q, uint16 **Kq);

// --------------------------- //
// --- fonctions de calcul --- //
// --------------------------- //

void convolution5x5_decimation_U16_I64_U16(uint16 **X, int height, int width, uint16 **K, int q, uint16 **Y);
void resample_I32(sint32 **X, int height, int width, sint32 **Y);


void makeBorder_I16(uint16 **T, int i0, int i1, int j0, int j1, int border);
void extendDataD2_I16(uint16 **X, int h, int w);
#ifdef __cplusplus
}
#endif

#endif
