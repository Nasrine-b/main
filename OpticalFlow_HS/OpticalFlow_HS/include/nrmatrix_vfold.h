/* ------------------------ */
/* --- nrmatrix_border_vfold.h --- */
/* ------------------------ */

/*
 * NRC matrix with vertical folded border
 */

/*
 * Copyright (c) 2017 Lionel Lacassagne, LIP6, UPMC, all rights reserved
 */

#ifndef __NR_matrix_border_vfold_H__
#define __NR_matrix_border_vfold_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

// --------------------------------------------- //
// -- internal functions: only for expert use -- //
// --------------------------------------------- //

uint8**   ui8matrix2_border_vfold(int height, int width, int border);
sint16** si16matrix2_border_vfold(int height, int width, int border);

sint32** si32matrix2_border_vfold(int height, int width, int border);
sint64** si64matrix2_border_vfold(int height, int width, int border);
float32** f32matrix2_border_vfold(int height, int width, int border);
float64** f64matrix2_border_vfold(int height, int width, int border);

uint8**   ui8matrix_border_vfold(int i0, int i1, int j0, int j1, int border);
sint16** si16matrix_border_vfold(int i0, int i1, int j0, int j1, int border);
sint32** si32matrix_border_vfold(int i0, int i1, int j0, int j1, int border);
sint64** si64matrix_border_vfold(int i0, int i1, int j0, int j1, int border);
float32** f32matrix_border_vfold(int i0, int i1, int j0, int j1, int border);
float64** f64matrix_border_vfold(int i0, int i1, int j0, int j1, int border);

void free_ui8matrix2_border_vfold (  uint8** m, int height, int width, int border);
void free_si16matrix2_border_vfold( sint16** m, int height, int width, int border);
void free_si32matrix2_border_vfold( sint32** m, int height, int width, int border);
void free_si64matrix2_border_vfold( sint64** m, int height, int width, int border);
void free_f32matrix2_border_vfold (float32** m, int height, int width, int border);
void free_f64matrix2_border_vfold (float64** m, int height, int width, int border);

void free_ui8matrix_border_vfold  ( uint8** m, int i0, int i1, int j0, int j1, int border);
void free_si16matrix_border_vfold (sint16** m, int i0, int i1, int j0, int j1, int border);
void free_si32matrix_border_vfold (sint32** m, int i0, int i1, int j0, int j1, int border);
void free_si64matrix_border_vfold (sint64** m, int i0, int i1, int j0, int j1, int border);
void free_f32matrix_border_vfold (float32** m, int i0, int i1, int j0, int j1, int border);
void free_f64matrix_border_vfold (float64** m, int i0, int i1, int j0, int j1, int border);

float32** f32matrix_vfold(int i0, int i1, int j0, int j1, int border);
void free_f32matrix_vfold (float32** m, int i0, int i1, int j0, int j1, int border);


void convert_ui8matrix_f32matrix_border_vfold(uint8** X, int i0, int i1, int j0, int j1, int border, float32 **Y);

void dup_f32matrix_border_vfold(float32** X, int i0, int i1, int j0, int j1, int border, float32 **Y);

#ifdef __cplusplus
}
#endif

#endif
