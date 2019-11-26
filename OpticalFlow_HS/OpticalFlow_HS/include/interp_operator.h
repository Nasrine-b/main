/* ------------------------- */
/* --- interp_operator.h --- */
/* ------------------------- */

/*
 * basic interpolation operator
 */

/*
 * Copyright (c) 2006-2017 Lionel Lacassagne, all rights reserved, LIP6, UPMC, CNRS
 */

#ifndef __INTERP_OPERATOR_H__
#define __INTERP_OPERATOR_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

float32 linear_interpolation(float32 f1, float32 f2, float32 x);

float32 cubic_interpolation_poly     (float32 f1, float32 f2, float32 f3, float32 f4, float32 x);
float32 cubic_interpolation_horner   (float32 f1, float32 f2, float32 f3, float32 f4, float32 x);
float32 cubic_interpolation_cell     (float32 f1, float32 f2, float32 f3, float32 f4, float32 x);

float32 cubic_interpolation_ipol     (float32 f1, float32 f2, float32 f3, float32 f4, float32 x);
float32 cubic_interpolation_burke    (float32 f0, float32 f1, float32 f2, float32 f3, float32 x);
float32 cubic_interpolation_breeuwsma(float32 f1, float32 f2, float32 f3, float32 f4, float32 x);

#ifdef __cplusplus
}
#endif

#endif // __INTERP_OPERATOR_H__

