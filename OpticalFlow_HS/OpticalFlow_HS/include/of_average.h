/* -------------------- */
/* --- OF_average.h --- */
/* -------------------- */

/*
 * Lukas Kanade main header
 */

/*
 * Copyright (c) 2006-2009 Lionel Lacassagne, IEF, all rights reserved
 */

#ifndef __OF_AVERAGE_H__
#define __OF_AVERAGE_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

// -----------------------------
// --- fonction de moyennage ---
// -----------------------------

// fonctions d'interface
void moyenneur_F32(float32** restrict X, int height, int width, int rayon, float32** restrict S, float32** restrict Y);
void moyenneur_F64(float32** restrict X, int height, int width, int rayon, float64** restrict S, float32** restrict Y);

// version filtres 2D: precis
void average2_F32(float32 **X, int height, int width, int rayon, float32 **Y);
void average2_F64(float32 **X, int height, int width, int rayon, float32 **Y);

// version filtres 1D: rapide
void average1_F32(float32 **X, int h, int w, int r, float32 **S, float32 **Y);


#ifdef __cplusplus
}
#endif

#endif
