#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "base_matrice.h"
#include "ope.h"
#include "mat_co_var.h"

/* fonction qui réalise un produit matricielle
*param A une matrice (pointeur)
*param B une autre matrice (pointeur)
*param nombre de ligne de A
*param nombre de ligne de B
*param nombre de colonne de A
*param nombre de colonne de B
*return un pointeur avec le produit des 2 autres matrices
*/
int **produit_M(int **A, int **B, int nb_l_A, int nb_l_B,int nb_c_A,int nb_c_B)
{
	//vérification faisabilité du produit AB
	if (nb_c_A != nb_l_B)
	{
		printf("Erreur de dimension\n");
		exit(1);
	}

	int **res;
	res = (int **) malloc( nb_l_A* sizeof(int *) );
	//vérification allocation mémoire
	if (res == NULL)
	{
		printf("Erreur allocation mémoire\n");
		exit(1);
	}
	for(int i=0; i<nb_l_A;i++)
	{
		res[i]=(int *)malloc(nb_c_B* sizeof(int));
		//vérification allocation mémoire
		if (res[i] == NULL)
		{
			printf("Erreur allocation mémoire\n");
			exit(1);
		}
	}


	
	//produit matricielle res=AB
	for (int i = 0; i < nb_l_A; i++)
	{
		for (int j = 0; j < nb_c_B; j++)
		{
			res[i][j]=0;
			for (int k = 0; k < nb_c_B; k++)
			{
				res[i][j]=res[i][j]+A[i][k]*B[k][j];
			}
		}
	}
	return res;
}

int main()
{
	//Déclaration des variables
	int dim=4;
	int ****vect=allocation_sous_mat(dim, dim);
	sous_mat_zeros(vect,dim,dim);
	affichage_sous_mat(vect,dim,dim);

	liberation_sous_mat(vect,dim,dim);

	return 0;
}