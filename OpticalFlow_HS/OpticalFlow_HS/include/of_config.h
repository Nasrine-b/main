/* ------------------- */
/* --- OF_config.h --- */
/* ------------------- */

/*
 * Lukas Kanade main header
 */

/*
 * Copyright (c) 2017 Lionel Lacassagne, all rights reserved
 * LIP6, UPMC, CNRS
 */

#ifndef __OF_CONFIG_H__
#define __OF_CONFIG_H__

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

#ifdef WIN32
#define FREQ 2.8e9
#else
#define FREQ 2.3e9
#endif

/*
 * activer ou desactiver les options suivantes
 * -------------------------------------------
 */

//#define OPENMP

// -------
#ifdef CLI
// -------
//#define ENABLE_OF_CHRONO
#define ENABLE_FAST
#define ENABLE_FAST2
    // a mettre dans le gestionnaire de projet
#define ENABLE_X

// pour le debug 2017 de HS
//#define SAVE_MODE
//#define SAVE_UV_MODE
#define SAVE_REC_MODE
#define SAVE_REC_FINAL_MODE
    
#define VERBOSE_MODE

#endif


// ========
#ifndef CLI
// ========

#define ENABLE_FAST
#define ENABLE_FAST2
// a mettre dans le gestionnaire de projet
//#define ENABLE_X
    
//#define ENABLE_OF_CHRONO
    
#define SAVE_MODE
#define SAVE_REC_MODE
#define SAVE_REC_FINAL_MODE

#define VERBOSE_MODE
    
//#define FINAL_DISPLAY_MODE
#define EXTEND_MODE
#define SAVE_UV_MODE
//#define SAVE_UV_FINAL_MODE
//#define ENABLE_F16
//#define ENABLE_F16_NORM
#define ENABLE_XP
    
#endif // !CLI
/*
 * ENABLE_FAST        : version rapide mais instable
 * ENABLE_X           : versions d'eXperimentations en cours
 
 * ENABLE_OF_CHRONO   : chronometrie
 * SAVE_MODE          : sauvegarde des donnees produites (image et matrice)
 * SAVE_REC_MODE      : sauvegarde des donnees reconstruites
 * SAVE_UV_MOD        : sauvegarde (txt et dat) des cartes (u,v) durant les iterations
 * SAVE_UV_FINAL_MODE : sauvegarde (txt et dat) de la carte (u,v) finale a la fin des iterations
 * VERBOSE_MODE       : mode verbeux
 * FINAL_DISPLAY_MODE : affichage uniquement a la fin
 * EXTEND_MODE        : extension des matrices (cf bug Matlab dans interpolation)
 * ENABLE_F16         : calcul en mode F16
 * ENABLE_F16_NORM    : troncature pour comparaison fairplay F32 et F16
 * ENABLE_XP          : eXtreme Programming: inclusion des fonctions de test
 */
    
// --------------------- //
// -- selecteurs SIMD -- //
// --------------------- //

// decommenter les ENABLE_SIMD au fur et a mesure que les fonctions sont recodees en SSE
// fait a l'epoque d'Altivec (AV)
//#define OF_ENABLE_SIMD_LKMTX
//#define OF_ENABLE_SIMD_CALC_IT
//#define OF_ENABLE_SIMD_CALC_GH
//#define OF_ENABLE_SIMD_CALC_UV
//#define OF_ENABLE_SIMD_INTERP
//#define OF_ENABLE_SIMD_GRAD
//#define OF_ENABLE_SIMD_ALLOC

// il faut au moins SIMD_ALLOC pour que les pointeurs soient des pointeurs SIMD
    
//#define OF_ENABLE_SIMD // pour pointeur / wrapper
    
/* ================================================ */
/* === ne rien ecrire en dessous de cette ligne === */
/* ================================================ */

#ifdef ENABLE_FAST
#define FAST(X) X
#else
#define FAST(X)
#endif

#if defined (OF_ENABLE_SIMD_ALLOC) || defined(OF_ENABLE_SIMD_LKMTX) || defined(OF_ENABLE_SIMD_CALC_IT) || defined(OF_ENABLE_SIMD_CALC_GH) || defined(OF_ENABLE_SIMD_CALC_UV) || defined(OF_ENABLE_SIMD_INTERP)
#define OF_ENABLE_SIMD
#endif
    
#ifdef ENABLE_OF_CHRONO
#ifdef WIN32
#pragma message("+++ OF Chrono enable +++")
#endif
#define OF_CHRONO(X) X
#define OF_CHRONO_VAR double t1,t2,dt
#define OF_CHRONO_START t1=dtime();
#define OF_CHRONO_STOP  t2=dtime(); dt=t2-t1
#else
#ifdef WIN32
#pragma message("--- OF Chrono disable ---")
#endif
#define OF_CHRONO(X)
#define OF_CHRONO_VAR
#define OF_CHRONO_START
#define OF_CHRONO_STOP
#endif

#ifdef SAVE_MODE
#define SAVE(X) X
#else
#define SAVE(X)
#endif

#ifdef SAVE_REC_MODE
#define SAVEREC(X) X
#else
#define SAVEREC(X)
#endif
// -----------------------
#ifdef SAVE_REC_FINAL_MODE
#define SAVEREC_FINAL(X) X
#else
#define SAVEREC_FINAL(X)
#endif
    
// ----------------
#ifdef SAVE_UV_MODE
#define SAVE_UV(X) X
#else
#define SAVE_UV(X)
#endif

// -----------------------
#ifdef SAVE_UV_FINAL_MODE
#define SAVE_UV_FINAL(X) X
#else
#define SAVE_UV_FINAL(X)
#endif

// ---------------
#ifdef EXTEND_MODE
#define OF_EXTEND(X) X
#else
#define OF_EXTEND(X)
#endif
    
// ----------------
#ifdef VERBOSE_MODE
#define OF_VERBOSE(X) X

// -----------------------
#ifdef FINAL_DISPLAY_MODE
#define FINAL_DISPLAY(X) X
#define CURRENT_DISPLAY(X)
#else
#define FINAL_DISPLAY(X)
#define CURRENT_DISPLAY(X) X
#endif
    
#else

#define OF_VERBOSE(X)
#define CURRENT_DISPLAY(X)
#define FINAL_DISPLAY(X)
#endif
    
/*
 * 256 empeche toute convergence ... (0.044)
 * 128 fait stabiliser l'algo a une erreur + grande (0.028)
 *  64 est OK mais moins bon (F32 et F16)
 *  32 est OK
 *  16 est OK
 *   8 est OK
 *   4 est OK
 *   2 est OK mais largement moins bon...
 *   1 fait stabiliser l'algo a une erreur + grande => normalisation necessaire ...
 */
#ifdef ENABLE_F16_NORM
#define F16NORM 32.0f
//#define F16NORM 1.0f
#else
#define F16NORM 1.0f
#endif

#ifdef __cplusplus
}
#endif

#endif
