/* -------------------- */
/* --- interp_fix.h --- */
/* -------------------- */

/*
 * fonction d'interpolation 2D en virgule fixe
 */

/*
 * Copyright (c) 2018-2018 Lionel Lacassagne, all rights reserved, Sorbonne University, LIP6, CNRS
 */

/*
 * 2018-07-10: creation interpCubic_I16_I32
 * 2018-12-11: ajout interpLinear_I16_I32
 */
#ifndef __INTERP_FIX_H__
#define __INTERP_FIX_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

void interpLinear_I16_I32(uint16 **I, int qi, sint32 **U, sint32 **V, int qs, int height, int width, uint16 **Irec);

void interpCubic_I16_I32(uint16 **I, int qi, sint32 **U, sint32 **V, int qs, int height, int width, uint16 **Irec);
void interpCubic_I16_F32(uint16 **I, int qi, sint32 **U, sint32 **V, int qs, int height, int width, uint16 **Irec);

#ifdef __cplusplus
}
#endif

#endif // __INTERP_FIX_H__

