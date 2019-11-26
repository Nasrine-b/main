/* -------------------- */
/* --- nr_pyramid.c --- */
/* -------------------- */

/*
 * NRC pyramid
 */

/*
 * Copyright (c) 2006-2009 Lionel Lacassagne, IEF, all rights reserved
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "of_config.h"

#include "nrc.h"
#include "nrdef.h"
#include "nrtype.h"
#include "nralloc.h"
#include "nrarith.h"
#include "nrio.h"

#ifdef FOLKI_ENABLE_SIMD
#ifdef CLI
#include "vnrcAV.h"
#else
#ifdef __APPLE__
#include <vNRC/vnrc.h>
#endif
#endif
#endif

#include "nrpyramid.h"
#include "nrmatrix_vfold.h"

/* --------------------- */
/* --- Sparse matrix --- */
/* --------------------- */

// ---------------------------------------------------------------------
uint8 ***ui8pyramid_sparse(int level, int height, int width, int border)
// ---------------------------------------------------------------------
{
    uint8 ***p;

    p = (uint8***) malloc(level * sizeof(uint8**));
    if(!p) {
        puts("Error1 allocating ui8pyramid...");
        exit(-1);
    }

    // for(int k=0; k<level; k++) { p[k] = 0; } // debug

    for(int k=0; k<level; k++) {
        p[k] = ui8matrix0(   0-border, height-1+border, 0-border, width-1+border);
        //zero_ui8matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    return p;
}
// -----------------------------------------------------------------------
sint16 ***si16pyramid_sparse(int level, int height, int width, int border)
// -----------------------------------------------------------------------
{
    sint16 ***p;
    
    p = (sint16***) malloc(level * sizeof(sint16**));
    if(!p) {
        puts("Error1 allocating si16pyramid...");
        exit(-1);
    }
    
    // for(int k=0; k<level; k++) { p[k] = 0; } // debug
    
    for(int k=0; k<level; k++) {
        p[k] = si16matrix( 0-border, height-1+border, 0-border, width-1+border);
        //zero_ui8matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    return p;
}
// -----------------------------------------------------------------------
uint16 ***ui16pyramid_sparse(int level, int height, int width, int border)
// -----------------------------------------------------------------------
{
    uint16 ***p;
    
    p = (uint16***) malloc(level * sizeof(uint16**));
    if(!p) {
        puts("Error1 allocating ui16pyramid...");
        exit(-1);
    }
    
    for(int k=0; k<level; k++) {
        p[k] = ui16matrix( 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    return p;
}
// -----------------------------------------------------------------------
sint32 ***si32pyramid_sparse(int level, int height, int width, int border)
// -----------------------------------------------------------------------
{
    sint32 ***p;
    
    p = (sint32***) malloc(level * sizeof(sint32**));
    if(!p) {
        puts("Error1 allocating si32pyramid...");
        exit(-1);
    }
    
    for(int k=0; k<level; k++) {
        p[k] = si32matrix( 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    return p;
}
// -----------------------------------------------------------------------
sint64 ***si64pyramid_sparse(int level, int height, int width, int border)
// -----------------------------------------------------------------------
{
    sint64 ***p;
    
    p = (sint64***) malloc(level * sizeof(sint64**));
    if(!p) {
        puts("Error1 allocating si64pyramid...");
        exit(-1);
    }
    
    // for(int k=0; k<level; k++) { p[k] = 0; } // debug
    
    for(int k=0; k<level; k++) {
        p[k] = si64matrix( 0-border, height-1+border, 0-border, width-1+border);
        //zero_ui8matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    return p;
}
// -----------------------------------------------------------------------
float32 ***f32pyramid_sparse(int level, int height, int width, int border)
// -----------------------------------------------------------------------
{
    // sparse allocation
    float32 ***p;

    p = (float32***) malloc(level * sizeof(float32**));
    if(!p) {
        puts("Error1 allocating f32pyramid...");
        exit(-1);
    }

    // for(int k=0; k<level; k++) { p[k] = 0; } // debug

    for(int k=0; k<level; k++) {
        p[k] = f32matrix0(0-border, height-1+border, 0-border, width-1+border);
        //zero_f32matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    return p;
}
// ------------------------------------------------------------------------
float64 ***f64pyramid_sparse(int level, int height, int width, int border)
// ------------------------------------------------------------------------
{
    // sparse allocation
    
    float64 ***p;
    
    p = (float64***) malloc(level * sizeof(float64**));
    if(!p) {
        puts("Error1 allocating f64pyramid...");
        exit(-1);
    }
    
    // for(int k=0; k<level; k++) { p[k] = 0; } // debug
    
    for(int k=0; k<level; k++) {
        p[k] = f64matrix0(0-border, height-1+border, 0-border, width-1+border);
        //zero_f64matrix(p[niv], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    return p;
}
// ----------------------------------------------------------------------------------
void free_ui8pyramid_sparse(uint8 ***p, int level, int height, int width, int border)
// ----------------------------------------------------------------------------------
{
    for(int k=0; k<level; k++) {
        //zero_ui8matrix(p[k], 0-border, height-1+border, 0-border, width-1+border); // debug
        free_ui8matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    free(p);
}
// ------------------------------------------------------------------------------------
void free_si16pyramid_sparse(sint16 ***p, int level, int height, int width, int border)
// ------------------------------------------------------------------------------------
{
    for(int k=0; k<level; k++) {
        free_si16matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    free(p);
}
// ------------------------------------------------------------------------------------
void free_ui16pyramid_sparse(uint16 ***p, int level, int height, int width, int border)
// ------------------------------------------------------------------------------------
{
    for(int k=0; k<level; k++) {
        free_ui16matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    free(p);
}
// ------------------------------------------------------------------------------------
void free_si32pyramid_sparse(sint32 ***p, int level, int height, int width, int border)
// ------------------------------------------------------------------------------------
{
    for(int k=0; k<level; k++) {
        free_si32matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    free(p);
}
// ------------------------------------------------------------------------------------
void free_si64pyramid_sparse(sint64 ***p, int level, int height, int width, int border)
// ------------------------------------------------------------------------------------
{
    for(int k=0; k<level; k++) {
        free_si64matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    free(p);
}
// ------------------------------------------------------------------------------------
void free_f32pyramid_sparse(float32 ***p, int level, int height, int width, int border)
// ------------------------------------------------------------------------------------
{
    // sparse allocation

    for(int k=0; k<level; k++) {
        //zero_f32matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        free_f32matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    free(p);
}
// ------------------------------------------------------------------------------------
void free_f64pyramid_sparse(float64 ***p, int level, int height, int width, int border)
// ------------------------------------------------------------------------------------
{
    // sparse allocation
    
    for(int k=0; k<level; k++) {
        //zero_f64matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        free_f64matrix(p[k], 0-border, height-1+border, 0-border, width-1+border);
        height = height / 2;
        width  = width  / 2;
    }
    free(p);
}

/* -------------------- */
/* --- Dense matrix --- */
/* -------------------- */

// -------------------------------------------------------------------
float32 ***f32pyramid_zeroBorder_v3(int level, int height, int width)
// -------------------------------------------------------------------
{
    // dense allocation
    int i, h, w, s;

    int asize; // accumulated size
    int aheight; // accumulated height
    int ai; // accumulated i

    float32 ***p;
    float32 **pp;
    float32 *ppp;

    // total size computation
    aheight = 0;
    asize   = 0;
    
    for(int k = 0; k < level; k++) {
        h = height >> k;
        w = width  >> k;
        s = h * w;

        aheight += h;
        asize   += s;
    }

    // _pyramid pointer allocation
    p = (float32***) malloc(level * sizeof(float32**));
    if(!p) {
        puts("Error1 allocating f32pyramid...");
        exit(-1);
    }

    // matrix pointer allocation
    pp = (float32**) malloc(aheight * sizeof(float32*));
    if(!pp) {
        puts("Error2 allocating f32pyramid...");
        exit(-1);
    }

    // data allocation
    ppp = (float32*) malloc(asize *sizeof(float32));
    if(!pp) {
        puts("Error3 allocating f32pyramid...");
        exit(-1);
    }
    pp[0] = ppp;
    p[0] = pp;

    // wrapping p to pp
    for(int k = 1; k < level; k++) {
        h = height >> (k-1);
        p[k] = p[k-1] + h;
    }

    // wrapping pp to ppp (pointer to first line of each matrix
    for(int k = 1; k < level; k++) {
        h = height >> (k-1);
        w = width  >> (k-1);
        s = h * w;
        pp[k] = pp[k-1] + s;
    }

    // wrapping data
    ai = 0;
    for(int k = 0; k < level; k++) {
        h = height >> k;
        w = width  >> k;
        s = h * w;
        for(i=1; i<h; i++) {
            pp[ai+i] = pp[ai+i-1] + w;
        }
        ai += h;   
    }

    return p;
}
// ----------------------------------------------------------------
float32 ***f32pyramid_zeroBorder(int level, int height, int width)
// ----------------------------------------------------------------
{
    // dense allocation

    int i, h, w, s;

    int asize; // accumulated size
    int aheight; // accumulated height
    int ai; // accumulated i

    float32 ***p;
    //float32 **pp;
    //float32 *ppp;

    // total size computation
    aheight = 0;
    asize  = 0;
    for(int k = 0; k < level; k++) {
        h = height >> k;
        w = width  >> k;
        s = h * w;

        aheight += h;
        asize   += s;
    }

    // _pyramid pointer allocation
    p = (float32***) malloc(level * sizeof(float32**));
    if(!p) {
        puts("Error1 allocating f32pyramid...");
        exit(-1);
    }

    // matrix pointer allocation
    p[0] = (float32**) malloc(aheight * sizeof(float32*));
    if(!p[0]) {
        puts("Error2 allocating f32pyramid...");
        exit(-1);
    }

    // data allocation
    p[0][0] = (float32*) malloc(asize *sizeof(float32));
    if(!p[0][0]) {
        puts("Error3 allocating f32pyramid...");
        exit(-1);
    }

    // wrapping p to pp
    for(int k = 1; k < level; k++) {
        h = height >> (k-1);
        p[k] = p[k-1] + h;
    }

    // wrapping pp to ppp (pointer to first line of each matrix
    for(int k = 1; k < level; k++) {
        h = height >> (k-1);
        w = width  >> (k-1);
        s = h * w;
        p[k][0] = p[k-1][0] + s;
    }

    // wrapping data
    ai = 0;
    for(int k = 0; k < level; k++) {
        h = height >> k;
        w = width  >> k;
        for(i=1; i<h; i++) {
            p[0][ai+i] = p[0][ai+i-1] + w;
        }
        ai += h;   
    }

    return p;
}
// ---------------------------------------------------------------------
uint8 ***ui8pyramid_dense(int level, int height, int width, int border)
// ---------------------------------------------------------------------
{
    // as border is a positive number, the offset should be added not substracted as for cube...
    
    int i, h, w; // dimensions of usefull data at a given level (without border)
    int  hb, wb, sb; // dimensions *with* border
    
    int asize; // accumulated size
    int aheight; // accumulated height
    //int ai; // accumulated i

    uint8 ***p;

    // total size computation
    aheight = 0;
    asize  = 0;
    for(int k = 0; k < level; k++) {
        h = height >> k; hb = h + 2 * border;
        w = width  >> k; wb = w + 2 * border;
        sb = hb * wb;
        
        aheight += hb;
        asize   += sb;
    }

    // _pyramid pointer allocation
    p = (uint8***) malloc(level * sizeof(uint8**));
    if(!p) {
        puts("Error1 allocating ui8pyramid...");
        exit(-1);
    }

    // matrix pointer allocation
    p[0] = (uint8**) malloc(aheight * sizeof(uint8*));
    if(!p[0]) {
        puts("Error2 allocating ui8pyramid...");
        exit(-1);
    }
    p[0] += border;

    // data allocation
    p[0][-border] = (uint8*) malloc(asize *sizeof(uint8));
    if(!p[0][-border]) {
        puts("Error3 allocating ui8pyramid...");
        exit(-1);
    }
    p[0][-border] += border;

    // wrapping p to pp
    for(int k = 1; k < level; k++) {
        h = height >> (k-1); hb = h + 2 * border;
        p[k] = p[k-1] + hb;
    }

    // wrapping pp to ppp (pointer to first line of each matrix
    for(int k = 1; k < level; k++) {
        h = height >> (k-1); hb = h + 2 * border;
        w = width  >> (k-1); wb = w + 2 * border;
        sb = hb * wb;
        p[k][-border] = p[k-1][-border] + sb;
    }

    // wrapping data
    //ai = 0;
    for(int k = 0; k < level; k++) {
        h = height >> level;
        w = width  >> level; wb = w + 2 * border;
        for(i=-border+1; i<h+border; i++) {
            p[k][i] = p[k][i-1] + wb;
        }
    }
    return p;
}
// -----------------------------------------------------------------------
float32 ***f32pyramid_dense(int level, int height, int width, int border)
// -----------------------------------------------------------------------
{
    // as border is a positive number, the offset should be added not substracted as for cube...
    
    int i, h, w; // dimensions of usefull data at a given level (without border)
    int  hb, wb, sb; // dimensions *with* border

    int asize; // accumulated size
    int aheight; // accumulated height

    float32 ***p;
    
    // total size computation
    aheight = 0;
    asize  = 0;
    for(int k = 0; k < level; k++) {
        h = height >> level; hb = h + 2 * border;
        w = width  >> level; wb = w + 2 * border;
        sb = hb * wb;

        aheight += hb;
        asize   += sb;
    }

    // _pyramid pointer allocation
    p = (float32***) malloc(level * sizeof(float32**));
    if(!p) {
        puts("Error1 allocating f32pyramid...");
        exit(-1);
    }

    // matrix pointer allocation
    p[0] = (float32**) malloc(aheight * sizeof(float32*));
    if(!p[0]) {
        puts("Error2 allocating f32pyramid...");
        exit(-1);
    }
    p[0] += border;

    // data allocation
    p[0][-border] = (float32*) malloc(asize *sizeof(float32));
    if(!p[0][-border]) {
        puts("Error3 allocating f32pyramid...");
        exit(-1);
    }
    p[0][-border] += border;

    // wrapping p to pp
    for(int k = 1; k < level; k++) {
        h = height >> (k-1); hb = h + 2 * border;
        p[k] = p[k-1] + hb;
    }

    // wrapping pp to ppp (pointer to first line of each matrix
    for(int k = 1; k < level; k++) {
        h = height >> (k-1); hb = h + 2 * border;
        w = width  >> (k-1); wb = w + 2 * border;
        sb = hb * wb;
        p[k][-border] = p[k-1][-border] + sb;
    }

    // wrapping data
    for(int k = 0; k < level; k++) {
        h = height >> k;
        w = width  >> k;
        wb = w + 2 * border;
        for(i=-border+1; i<h+border; i++) {
            p[k][i] = p[k][i-1] + wb;
        }
    }
    return p;
}
// ----------------------------------------------------------------------------------
void free_ui8pyramid_dense(uint8 ***p, int level, int height, int width, int border)
// ----------------------------------------------------------------------------------
{
    puts("free_ui8pyramid_dense NOT YET IMPLEMENTED");
    exit(-1);
}
// ------------------------------------------------------------------------------------
void free_f32pyramid_dense(float32 ***p, int level, int height, int width, int border)
// ------------------------------------------------------------------------------------
{
    puts("free_f32pyramid_dense NOT YET IMPLEMENTED");
    exit(-1);
}
// ------------------------------------------------------------------------------------
void free_f64pyramid_dense(float64 ***p, int level, int height, int width, int border)
// ------------------------------------------------------------------------------------
{
    puts("free_f64pyramid_dense NOT YET IMPLEMENTED");
    exit(-1);
}

/* ------------------------ */
/* --- Public functions --- */
/* ------------------------ */

// --------------------------------------------------------------
uint8 ***ui8pyramid(int level, int height, int width, int border)
// --------------------------------------------------------------
{
#ifdef DENSE_PYRAMID_ALLOCATION
    return ui8pyramid_dense(level, height, width, border);
#else
    return ui8pyramid_sparse(level, height, width, border);
#endif
}
// ----------------------------------------------------------------
sint16 ***si16pyramid(int level, int height, int width, int border)
// ----------------------------------------------------------------
{
#ifdef DENSE_PYRAMID_ALLOCATION
    return si16pyramid_dense(level, height, width, border);
#else
    return si16pyramid_sparse(level, height, width, border);
#endif
}
// ----------------------------------------------------------------
uint16 ***ui16pyramid(int level, int height, int width, int border)
// ----------------------------------------------------------------
{
#ifdef DENSE_PYRAMID_ALLOCATION
    return ui16pyramid_dense(level, height, width, border);
#else
    return ui16pyramid_sparse(level, height, width, border);
#endif
}
// ----------------------------------------------------------------
sint32 ***si32pyramid(int level, int height, int width, int border)
// ----------------------------------------------------------------
{
#ifdef DENSE_PYRAMID_ALLOCATION
    return si32pyramid_dense(level, height, width, border);
#else
    return si32pyramid_sparse(level, height, width, border);
#endif
}

// ----------------------------------------------------------------
sint64 ***si64pyramid(int level, int height, int width, int border)
// ----------------------------------------------------------------
{
#ifdef DENSE_PYRAMID_ALLOCATION
    return si64pyramid_dense(level, height, width, border);
#else
    return si64pyramid_sparse(level, height, width, border);
#endif
}
// ----------------------------------------------------------------
float32 ***f32pyramid(int level, int height, int width, int border)
// ----------------------------------------------------------------
{  
#ifdef DENSE_PYRAMID_ALLOCATION
    return f32pyramid_dense(level, height, width, border);
#else
    return f32pyramid_sparse(level, height, width, border);
#endif
}
// -----------------------------------------------------------------
float64 ***f64pyramid(int level, int height, int width, int border)
// -----------------------------------------------------------------
{  
#ifdef DENSE_PYRAMID_ALLOCATION
    return f64pyramid_dense(level, height, width, border);
#else
    return f64pyramid_sparse(level, height, width, border);
#endif
}
// ---------------------------------------------------------------------------
void free_ui8pyramid(uint8 ***p, int level, int height, int width, int border)
// ---------------------------------------------------------------------------
{

#ifdef DENSE_PYRAMID_ALLOCATION
    free_ui8pyramid_dense(p, level, height, width, border);
#else
    free_ui8pyramid_sparse(p, level, height, width, border);
#endif
}
// -----------------------------------------------------------------------------
void free_si16pyramid(sint16 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    
#ifdef DENSE_PYRAMID_ALLOCATION
    free_si16pyramid_dense(p, level, height, width, border);
#else
    free_si16pyramid_sparse(p, level, height, width, border);
#endif
}
// -----------------------------------------------------------------------------
void free_ui16pyramid(uint16 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    
#ifdef DENSE_PYRAMID_ALLOCATION
    free_ui16pyramid_dense(p, level, height, width, border);
#else
    free_ui16pyramid_sparse(p, level, height, width, border);
#endif
}
// -----------------------------------------------------------------------------
void free_si32pyramid(sint32 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    
#ifdef DENSE_PYRAMID_ALLOCATION
    free_si32pyramid_dense(p, level, height, width, border);
#else
    free_si32pyramid_sparse(p, level, height, width, border);
#endif
}
// -----------------------------------------------------------------------------
void free_si64pyramid(sint64 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    
#ifdef DENSE_PYRAMID_ALLOCATION
    free_si64pyramid_dense(p, level, height, width, border);
#else
    free_si64pyramid_sparse(p, level, height, width, border);
#endif
}
// -----------------------------------------------------------------------------
void free_f32pyramid(float32 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
#ifdef DENSE_PYRAMID_ALLOCATION
    free_f32pyramid_dense(p, level, height, width, border);
#else
    free_f32pyramid_sparse(p, level, height, width, border);
#endif
}
// -----------------------------------------------------------------------------
void free_f64pyramid(float64 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
#ifdef DENSE_PYRAMID_ALLOCATION
    free_f64pyramid_dense(p, level, height, width, border);
#else
    free_f64pyramid_sparse(p, level, height, width, border);
#endif
}

/* ----------------- */
/* --- Functions --- */
/* ----------------- */

// ----------------------------------------------------------------------------
void zero_f32pyramid_zeroBorder(float32 ***p, int level, int height, int width)
// ----------------------------------------------------------------------------
{
    int h, w;
    int i, j;

    float32 **pp;
    float32 *ppp;

    pp  = p[0];
    ppp = pp[0];
    for(int k  =0; k < level; k++) {

        h = height >> k;
        w = width  >> k;
             
        for(i=0; i<h; i++) {
            for(j=0; j<w; j++) {
                *ppp++ = 0.0f;
            }
        }
    }
}
// ---------------------------------------------------------------------------
void zero_ui8pyramid(uint8 ***p, int level, int height, int width, int border)
// ---------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        
        int h = height >> k;
        int w = width  >> k;
        
        for(int i=0-border; i<h+border; i++) {
            for(int j=0-border; j<w+border; j++) {
                p[k][i][j] = 0;
            }
        }
    }
}
// -----------------------------------------------------------------------------
void zero_si16pyramid(sint16 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        
        int h = height >> k;
        int w = width  >> k;
        
        for(int i=0-border; i<h+border; i++) {
            for(int j=0-border; j<w+border; j++) {
                p[k][i][j] = 0;
            }
        }
    }
}
// -----------------------------------------------------------------------------
void zero_ui16pyramid(uint16 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        
        int h = height >> k;
        int w = width  >> k;
        
        for(int i=0-border; i<h+border; i++) {
            for(int j=0-border; j<w+border; j++) {
                p[k][i][j] = 0;
            }
        }
    }
}
// -----------------------------------------------------------------------------
void zero_si32pyramid(sint32 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        
        int h = height >> k;
        int w = width  >> k;
        
        for(int i=0-border; i<h+border; i++) {
            for(int j=0-border; j<w+border; j++) {
                p[k][i][j] = 0;
            }
        }
    }
}
// -----------------------------------------------------------------------------
void zero_si64pyramid(sint64 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        
        int h = height >> k;
        int w = width  >> k;
        
        for(int i=0-border; i<h+border; i++) {
            for(int j=0-border; j<w+border; j++) {
                p[k][i][j] = 0;
            }
        }
    }
}
// -----------------------------------------------------------------------------
void zero_f32pyramid(float32 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {

        int h = height >> level;
        int w = width  >> level;
             
        for(int i = 0-border; i < h+border; i++) {
            for(int j = 0-border; j < w+border; j++) {
                p[k][i][j] = 0.0f;
            }
        }
    }
}
// -----------------------------------------------------------------------------
void zero_f64pyramid(float64 ***p, int level, int height, int width, int border)
// -----------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        
        int h = height >> k;
        int w = width  >> k;
        
        for(int i = 0-border; i < h+border; i++) {
            for(int j = 0-border; j < w+border; j++) {
                p[k][i][j] = 0.0;
            }
        }
    }
}
// ---------------------------------------------------------------------------------------
void dup_ui8pyramid(uint8 ***X, int level, int height, int width, int border, uint8 ***Y)
// ---------------------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        int h = height >> k;
        int w = width  >> k;
        dup_ui8matrix(X[k], 0-border, h-1+border, 0-border, w-1+border, Y[k]);
    }
}
// -----------------------------------------------------------------------------------------
void dup_si16pyramid(sint16 ***X, int level, int height, int width, int border, sint16 ***Y)
// -----------------------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        int h = height >> k;
        int w = width  >> k;
        dup_si16matrix(X[k], 0-border, h-1+border, 0-border, w-1+border, Y[k]);
    }
}
// -----------------------------------------------------------------------------------------
void dup_ui16pyramid(uint16 ***X, int level, int height, int width, int border, uint16 ***Y)
// -----------------------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        int h = height >> k;
        int w = width  >> k;
        dup_ui16matrix(X[k], 0-border, h-1+border, 0-border, w-1+border, Y[k]);
    }
}
// -----------------------------------------------------------------------------------------
void dup_si32pyramid(sint32 ***X, int level, int height, int width, int border, sint32 ***Y)
// -----------------------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        int h = height >> k;
        int w = width  >> k;
        dup_si32matrix(X[k], 0-border, h-1+border, 0-border, w-1+border, Y[k]);
    }
}
// -----------------------------------------------------------------------------------------
void dup_si64pyramid(sint64 ***X, int level, int height, int width, int border, sint64 ***Y)
// -----------------------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        int h = height >> k;
        int w = width  >> k;
        dup_si64matrix(X[k], 0-border, h-1+border, 0-border, w-1+border, Y[k]);
    }
}
// -------------------------------------------------------------------------------------------
void dup_f32pyramid(float32 ***X, int level, int height, int width, int border, float32 ***Y)
// -------------------------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        int h = height >> level;
        int w = width  >> level;
        dup_f32matrix(X[k], 0-border, h-1+border, 0-border, w-1+border, Y[k]);
    }
}
// -------------------------------------------------------------------------------------------
void dup_f64pyramid(float64 ***X, int level, int height, int width, int border, float64 ***Y)
// -------------------------------------------------------------------------------------------
{
    for(int k = 0; k < level; k++) {
        int h = height >> k;
        int w = width  >> k;
        dup_f64matrix(X[k], 0-border, h-1+border, 0-border, w-1+border, Y[k]);
    }
}
