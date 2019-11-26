/* ----------------------------- */
/* --- HornShunck_operator.c --- */
/* ----------------------------- */

/*
 * Copyright (c) 2017, Lionel Lacassagne, all rights reserved
 * LIP6, UPMC, CNRS
 */

#ifndef _HORNSHUNCK_OPERATOR_H_
#define _HORNSHUNCK_OPERATOR_H_

#ifdef __cplusplus
extern "C" {
#endif

// flags pour le calcul en virgule fixe
#define ENABLE_GRADIENT_DIV
#define ENABLE_AVERAGE_DIV

// -------------------------
// mono-resolution functions
// -------------------------

void hs_basic_F32(float32 **E0, float32 **E1, float32 **U0, float32 **V0, int height, int width, float32 alpha2, float32 **U1, float32 **V1);

void hs_calc_gradientX_F32  (float32 **E0, float32 **E1, int height, int width, float32 **Ex);
void hs_calc_gradientY_F32  (float32 **E0, float32 **E1, int height, int width, float32 **Ey);
void hs_calc_gradientT_F32  (float32 **E0, float32 **E1, int height, int width, float32 **Et);
void hs_calc_gradientXYT_F32(float32 **E0, float32 **E1, int height, int width, float32 **Ex, float32 **Ey, float32 **Et);

void hs_calc_UmVm_F32(float32 **U, float32 **V, int height, int width, float32 **Um, float32 **Vm);

void hs_calc_UV0_F32      (float32 **Ex, float32 **Ey, float32 **Et,                             int height, int width, float32 alpha2, float32 **U, float32 **V);
void hs_calc_UV_F32       (float32 **Ex, float32 **Ey, float32 **Et, float32 **Um, float32 **Vm, int height, int width, float32 alpha2, float32 **U, float32 **V);
void hs_calc_UV_direct_F32(float32 **Ex, float32 **Ey, float32 **Et,                             int height, int width, float32 alpha2, float32 **U, float32 **V);

// --------------------------------------
// multi-resolution / pyramidal functions
// --------------------------------------

void hs_accumulate_UV_F32(float32 **dU, float32 **dV, int height, int width, float32 **U, float32 **V);
void hs_add_UV_F32(float32 **U, float32 **V, float32 **dU, float32 **dV, int height, int width, float32 **Ut, float32 **Vt);

// ----------------
// special function
// ----------------

void hs_clip_UV_F32(float32 **U, float32 **V, int height, int width);

// -----------------------
// -- High-level functions
// -----------------------

void hs_step0_F32(float32 **E0, float32 **E1, int height, int width, float32 alpha2, float32 **Ex, float32 **Ey, float32 **Et, float32 **U1, float32 **V1);
void hs_1step_F32(float32 **U0, float32 **V0, float32 **Ex, float32 **Ey, float32 **Et, int height, int width, float32 alpha2, float32 **Um, float32 **Vm, float32 **U1, float32 **V1);

#ifdef __cplusplus
}
#endif

#endif // _HORNSHUNCK_OPERATOR_H_