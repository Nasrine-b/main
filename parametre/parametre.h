#ifndef __PARAMETRE_H
#define __PARAMETRE_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float f_S(float* vec,int taille,int p);

float f_C(float* vec,int taille,int p);

float f_T(float Cp_b, float Sp_b);

float f_R(float C, float S);

float f_V(float R);

float f_sigma(float R);

float f_delta(float T,float R);

float f_s(float R1, float R2, float T1, float T2);

float f_k(float R1, float R2, float T1, float T2);

void circ_data(float* vec, int taille, float* R1, float* R1_b, float* theta, float* V, float* sigma, float* delta, float* s, float* k); 

#endif