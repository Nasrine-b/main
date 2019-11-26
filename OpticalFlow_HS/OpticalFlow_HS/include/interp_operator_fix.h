/* ----------------------------- */
/* --- interp_operator_fix.h --- */
/* ----------------------------- */

/*
 * interpmolation cubique en virgule fixe
 */

/*
 * Copyright (c) 2018 Lionel Lacassagne, all rights reserved, LIP6, UPMC, CNRS
 * 2018-12-09: creation
 */

#ifndef __INTERP_OPERATOR_FIX_H__
#define __INTERP_OPERATOR_FIX_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

sint64 linear_interpolation_I64(sint64 e1, sint64 e2, sint64 X, int q);

sint64 cubic_interpolation_ipol_I64(sint64 e1, sint64 e2, sint64 e3, sint64 e4, sint64 X, int q);
float32 cubic_interpolation_ipol_I16_F32(float32 e1, float32 e2, float32 e3, float32 e4, float32 X, int q);


#ifdef __cplusplus
}
#endif

#endif // __INTERP_OPERATOR_FIX_H__

