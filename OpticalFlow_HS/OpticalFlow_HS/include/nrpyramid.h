/* ------------------- */
/* --- nrpyramid.h --- */
/* ------------------- */

/*
 * NRC pyramid
 */

/*
 * Copyright (c) 2006-2009 Lionel Lacassagne, IEF, all rights reserved
 */

/*
 * height and width: size of the full resolution image (or data)
 * nbLevel: number of pyramid's levels, from 0 to nbLevel-1 [0:nbLevel-1]
 */

#ifndef __NR_PYRAMID_H__
#define __NR_PYRAMID_H__

//#define DENSE_PYRAMID_ALLOCATION

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

// --------------------------------------------- //
// -- internal functions: only for expert use -- //
// --------------------------------------------- //
     
uint8    ***ui8pyramid_sparse(int level, int height, int width, int border);
sint16  ***si16pyramid_sparse(int level, int height, int width, int border);
uint16  ***ui16pyramid_sparse(int level, int height, int width, int border);
sint32  ***si32pyramid_sparse(int level, int height, int width, int border);
float32  ***f32pyramid_sparse(int level, int height, int width, int border);
float64  ***f64pyramid_sparse(int level, int height, int width, int border);

void free_ui8pyramid_sparse (uint8   ***p, int level, int height, int width, int border);
void free_si16pyramid_sparse (sint16 ***p, int level, int height, int width, int border);
void free_ui16pyramid_sparse (uint16 ***p, int level, int height, int width, int border);
void free_si32pyramid_sparse (sint32 ***p, int level, int height, int width, int border);
void free_f32pyramid_sparse (float32 ***p, int level, int height, int width, int border);
void free_f64pyramid_sparse (float64 ***p, int level, int height, int width, int border);

uint8   ***ui8pyramid_dense(              int level, int height, int width, int border);
float32 ***f32pyramid_dense(              int level, int height, int width, int border);
void free_ui8pyramid_dense (uint8   ***p, int level, int height, int width, int border);
void free_f32pyramid_dense (float32 ***p, int level, int height, int width, int border);

void zero_f32pyramid_zeroBorder(float32 ***p, int level, int height, int width);

void zero_ui8pyramid_dense(uint8   ***p, int level, int height, int width, int border);
void zero_f32pyramid_dense(float32 ***p, int level, int height, int width, int border);
    
// ---------------------- //
// -- public functions -- //
// ---------------------- //

 uint8  *** ui8pyramid(int level, int height, int width, int border);
 sint16 ***si16pyramid(int level, int height, int width, int border);
 uint16 ***ui16pyramid(int level, int height, int width, int border);
 sint32 ***si32pyramid(int level, int height, int width, int border);
 sint64 ***si64pyramid(int level, int height, int width, int border);
float32 *** f32pyramid(int level, int height, int width, int border);
float64 *** f64pyramid(int level, int height, int width, int border);

void zero_ui8pyramid (uint8   ***p, int level, int height, int width, int border);
void zero_si16pyramid(sint16  ***p, int level, int height, int width, int border);
void zero_ui16pyramid(uint16  ***p, int level, int height, int width, int border);
void zero_si32pyramid(sint32  ***p, int level, int height, int width, int border);
void zero_si64pyramid(sint64  ***p, int level, int height, int width, int border);
void zero_f32pyramid (float32 ***p, int level, int height, int width, int border);
void zero_f64pyramid (float64 ***p, int level, int height, int width, int border);

void free_ui8pyramid (uint8   ***p, int level, int height, int width, int border);
void free_si16pyramid(sint16  ***p, int level, int height, int width, int border);
void free_ui16pyramid(uint16  ***p, int level, int height, int width, int border);
void free_si32pyramid(sint32  ***p, int level, int height, int width, int border);
void free_si64pyramid(sint64  ***p, int level, int height, int width, int border);
void free_f32pyramid (float32 ***p, int level, int height, int width, int border);
void free_f64pyramid (float64 ***p, int level, int height, int width, int border);

void display_ui8pyramid (uint8   ***p, int level, int height, int width, int border, char *format, char *name);
void display_si16pyramid(sint16  ***p, int level, int height, int width, int border, char *format, char *name);
void display_si64pyramid(sint64  ***p, int level, int height, int width, int border, char *format, char *name);
void display_f32pyramid (float32 ***p, int level, int height, int width, int border, char *format, char *name);
void display_f64pyramid (float64 ***p, int level, int height, int width, int border, char *format, char *name);

void dup_ui8pyramid (uint8   ***X, int level, int height, int width, int border, uint8   ***Y);
void dup_si16pyramid(sint16  ***X, int level, int height, int width, int border, sint16  ***Y);
void dup_ui16pyramid(uint16  ***X, int level, int height, int width, int border, uint16  ***Y);
void dup_si64pyramid(sint64  ***X, int level, int height, int width, int border, sint64  ***Y);
void dup_f32pyramid (float32 ***X, int level, int height, int width, int border, float32 ***Y);
void dup_f64pyramid (float64 ***X, int level, int height, int width, int border, float64 ***Y);

#ifdef __cplusplus
}
#endif

#endif
