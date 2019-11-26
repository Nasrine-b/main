/* -------------- */
/* --- main.c --- */
/* -------------- */

/*
 * Copyright (c) 2006-2009 Lionel Lacassagne, IEF, all rights reserved
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

// NRC lib
#include "nrc.h"
#include "nrdef.h"
#include "nrtype.h"
#include "nralloc.h"
#include "nrarith.h"
#include "nrio.h"
//#include "sequence.h"

//#include "of_config.h"
//#include "of.h"
//#include "folki_iter.h"
//#include "hs_iter.h"

//#include "bench_hs.h"
//#include "bench_folki.h"

#include "test_hs_operator.h"
//#include "test_lk_operator.h"
//#include "test_folki_operator.h"

//#include "bench_hs_operator.h"

#ifdef ENABLE_FAST
//#include "test_folki_function.h"
//#include "test_folki_interp.h"
#endif

/* ------------------- */
void display_DEFINE(void)
/* ------------------- */
{
    puts("#########################################");
    puts("### configuration #######################");
    puts("#########################################");
    
#ifdef ENABLE_FOLKI_CHRONO
    puts("ENABLE_FOLKI_CHRONO");
#endif
    
#ifdef SAVE_MODE
    puts("SAVE_MODE on");
#endif
    
#ifdef SAVE_REC_MODE
    puts("SAVE_REC_MODE on");
#endif
    
#ifdef VERBOSE_MODE
    puts("VERBOSE_MODE on");
#endif

#ifdef FINAL_DISPLAY_MODE
    puts("FINAL_DISPLAY_MODE on");
#endif
    
#ifdef EXTEND_MODE
    puts("EXTEND_MODE on");
#endif
    
#ifdef ENABLE_F16_NORM
    puts("ENABLE_F16_NORM");
#endif
 
#ifdef FOLKI_ENABLE_SIMD_LKMTX
    puts("FOLKI_ENABLE_SIMD_LKMTX");
#endif
#ifdef FOLKI_ENABLE_SIMD_CALC_IT
    puts("FOLKI_ENABLE_SIMD_CALC_IT");
#endif
#ifdef FOLKI_ENABLE_SIMD_CALC_GH
    puts("FOLKI_ENABLE_SIMD_CALC_GH");
#endif
#ifdef FOLKI_ENABLE_SIMD_CALC_UV
    puts("FOLKI_ENABLE_SIMD_CALC_UV");
#endif
    
#ifdef FOLKI_ENABLE_SIMD_INTERP
    puts("FOLKI_ENABLE_SIMD_INTERP");
#endif

#ifdef FOLKI_ENABLE_SIMD       
    puts("FOLKI_ENABLE_SIMD");
#endif

#ifdef FOLKI_ENABLE_FOLKI
    puts("FOLKI_ENABLE_FOLKI");
#endif

#ifdef ENABLE_X
    puts("ENABLE_X");
#endif

#ifdef ENABLE_FAST
    puts("ENABLE_FAST");
#endif

#ifdef ENABLE_FAST2
    puts("ENABLE_FAST2");
#endif
    
#ifdef CLI
    puts("CLI");
#endif
    
}
/* -------------------------------- */
int mainScript(int argc, char *argv[])
/* -------------------------------- */
{
    int i;
    
    /* ------------- */
    /* --- batch --- */
    /* ------------- */
    //display_DEFINE();
    
    for(i=0; i<argc; i++) {
        printf("[main]: %s \n", argv[i]);
    } puts("");  
    
#ifdef OF_ENABLE_FOLKI
    puts("[Main]: launching FOLKI_bench_sequence");
    FOLKI_bench_sequence(argc, argv);
#endif
    
#ifdef OF_ENABLE_HS
    puts("[Main]: launching HS_bench_sequence");
    HS_bench_sequence(argc, argv);
#endif
    return 0;
}
/* ------------------------------ */
int main_CLI(int argc, char *argv[])
/* ------------------------------ */
{
    puts("main_CLI");
    //FOLKI_test(); return 0; // test unitaire de LK sur 2 images
    
    //HS_test(); //return 0; // test de hs_iter
    
    HS_test_operator(argc, argv); return 0; // <- celui la
    //LK_test_operator(argc, argv); return 0; // <- celui la
    
    //LK_bench_withAlloc(); return 0;
    
    //main_bench_hs_operator(); return 0; // 2018
    
    //return mainScript(argc, argv);
}
/* ------------------------------------- */
int main_interactif(int argc, char *argv[])
/* ------------------------------------- */
{
    puts("main_interactif");
    
    //puts("il manque certaines fonctions HS_xxx_F16 qui ne sont pas codees");
    //puts("hello world\n");
    //for(i=0; i<argc; i++) printf("%s ", argv[i]); puts("\n");
    //puts("main: go");
    
#ifdef ENABLE_XP    
    //FOLKI_test();
    //HS_test();
    //HS_test_Fast(); return;
    
    //HS_test_operator(argc, argv);
    
     HS_test_operator(argc, argv);
     //LK_test_operator(argc, argv);
     //FOLKI_test_operator(argc, argv);
    
    //FOLKI_testTroncature();

    //testInterp();
    //testPrecision();
    //testPrecision2();
    //testMoyenneurKxK();
    //testCenteredSumKxK(); // somme cumulee, version speciale Onera
    //benchAverage(); // benchmark de vitesse
    //testScan(); // version avec somme cumulee sigma sliding
    //testAverage1();

    //testF16();
    //FOLKI_test();
    

    //testMultiple();
    //testPyramid();
    //test_vPyramid();
    //testModulo();
    //testAngle(); // pour le calcul de l'angle (second critere Onera)
    return 0;
    //return -2;
#endif

    // ----------------- //
    // --- benchmark --- //
    // ----------------- //

    //FOLKI_bench_withAlloc(); // probleme residuel de malloc ?
    //FOLKI_bench_withoutAlloc(); // <- this one
    //HS_bench_withoutAlloc();
    //HS_bench_withAlloc(); // <- OK pour Horn & Shunk
    
    // ------------------------ //
    // --- batch processing --- //
    // ------------------------ //
    
    //display_DEFINE();
    //mainScript(argc, argv);
    
    //FOLKI_bench_sequence(argc, argv); // <- this one
    //HS_bench_sequence(argc, argv); // <- this one
    
    //puts("main: bye bye");
    return 0;
}
/* ========================== */
int main(int argc, char *argv[])
/* ========================== */
{
    display_DEFINE();

#ifdef WIN32
    // en attendant de faire propre car il faut definir CLI pour que les headers soient vus.
    //testInterp(); return 0; // pour tester SSE4 en 2010
    //main_CLI(argc, argv);
    //FOLKI_bench_sequence(argc, argv);
    return 0;
#endif

#ifdef CLI
    return main_CLI(argc, argv);
#else
    return main_interactif(argc, argv);
#endif
}