/* ---------------- */
/* --- interp.h --- */
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

#ifndef __INTERP_H__
#define __INTERP_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

/* --------------------- */
/* --- interpolation --- */
/* --------------------- */

void interpNearest_F32(float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec);
void interpLinear_F32 (float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec);
void interpCubic_F32  (float32 **I, float32 **U, float32 **V, int height, int width, float32 **Irec);

void interpCubic_I64(uint8 **I, sint64 **U, sint64 **V, int height, int width, int q, uint8 **Irec);

/* ---------------- */
/* --- gradient --- */
/* ---------------- */

void gradientCubic_F32(float32 **I, float32 **U, float32 **V, int height, int width, float32 **Ix, float32 **Iy);

#ifdef __cplusplus
}
#endif

#endif // __INTERP_H__

