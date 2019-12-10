#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "base_matrice.h"

// void gen_cercle(norme, int sigma, int nb_vect)
// {
//   **res = mat_zeros(nb_vect, 2);
//   for(int i=0; i<nb_vect; i++) 
//   {

//     res[i](i,:) = [normrnd(norme, sigma), unifrnd(0, 2*pi())];
//   }
//   return res;
  
// }

int unif(int borne_inf, int borne_sup)
{
	return rand()%(borne_sup-borne_inf)+borne_inf;
}

int **runif(int borne_inf, int borne_sup, int taille)
{
	int **m;

	m=allocation(taille,1);
	for (int i = 0; i < taille; ++i)
	{
		m[i][0]=rand()%(borne_sup-borne_inf)+borne_inf;
	}
	return m;
}

int **rnorm(int mu, int sigma, int *x, int taille)
{

	int **m;

	m=allocation(taille,1);
	for (int i = 0; i < taille; ++i)
	{
		m[i][0]=1/(sigma*sqrt(2*M_PI)) * exp(-( pow(x[i]-mu,2)/(2*pow(sigma,2))) );
	}

	return m;
}

int power_int(int x, int p)
{
	if(p==0)
	{
		return 1;
	}
	if (p==1)
	{
		return x; 
	}
	return x*power_int(x,p-1);
}