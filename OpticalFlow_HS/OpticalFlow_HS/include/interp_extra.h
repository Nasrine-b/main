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

#ifndef __INTERP_EXTRA_H__
#define __INTERP_EXTRA_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

/* --------------------- */
/* --- interpolation --- */
/* --------------------- */
void initInterpCubicLUT_F32(float32 **lut, int32 size);
void initGradientCubicLUT_F32(float32 **lut, int32 size);

void interpLinear_I8 (uint8    **I, int height, int width, float32 **U, float32 **V, uint8    **Irec);
void interpLinear_I32(sint32  **I, int height, int width, float32 **U, float32 **V, sint32  **Irec);

void interpLacas_I8 (uint8    **I, int height, int width, float32 **U, float32 **V, uint8    **Irec);
void interpLacas_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **Irec);

void interpCubic_I8 (uint8    **I, int height, int width, float32 **U, float32 **V, uint8    **Irec);
void interpCubic_I32(sint32  **I, int height, int width, float32 **U, float32 **V, sint32  **Irec);

void interpCubic_F32_F64_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **Irec);

void interpCubicLUT_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **C, int32 size, float32 **Irec);

/* ---------------- */
/* --- gradient --- */
/* ---------------- */

void sobel_I8F32    (uint8   **I, int height, int width, float32 **Ix, float32 **Iy);
void sobel_F32     (float32 **I, int height, int width, float32 **Ix, float32 **Iy);
void gradient_I8F32(uint8    **I, int height, int width, float32 **Ix, float32 **Iy);

void gradientLinear_I8 (uint8    **I, int height, int width, float32 **U, float32 **V, float32 **Ix, float32 **Iy);
void gradientLinear_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **Ix, float32 **Iy);

void gradientCubic_I8 (uint8    **I, int height, int width, float32 **U, float32 **V, float32 **Ix, float32 **Iy);
void gradientCubic0_F32(float32 **I, int height, int width, float32 **Ix, float32 **Iy);
void gradientCubicLUT_F32(float32 **I, int height, int width, float32 **U, float32 **V, float32 **lutI, int32 sizeI, float32 **lutG, int32 sizeG, float32 **Ix, float32 **Iy);

void gradientLacas_I8(uint8 **I, int height, int width, float32 **U, float32 **V, float32 **Ix, float32 **Iy);

#ifdef __cplusplus
}
#endif

#endif // __INTERP_EXTRA_H__

