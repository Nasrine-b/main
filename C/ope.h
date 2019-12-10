#ifndef __ope_H_
#define __ope_H_

int unif(int borne_inf, int borne_sup);
int **runif(int borne_inf, int borne_sup, int taille);
int **rnorm(int mu, int sigma, int *x, int taille);
int power_int(int x, int p);

#endif