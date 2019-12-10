#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "base_matrice.h"
#include "ope.h"

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
	int **M1;
	int ligne_M1=3;
	int colonne_M1=3;
	int unif_test=unif(3,5);
	//------------------------

	//Affectation M1
	// printf("Nombre de ligne :");
	// scanf("%d",&ligne_M1);
	// while (ligne_M1 <= 0)
	// {
	// 	printf("Nombre de lignes incorrect\n Nombre de ligne:");
	// 	scanf("%d",&ligne_M1);
	// }

	// printf("Nombre de colonne :");
	// scanf("%d",&colonne_M1);
	// while (colonne_M1 <= 0)
	// {
	// 	printf("Nombre de colonne incorrect\n Nombre de colonne:");
	// 	scanf("%d",&colonne_M1);
	// }
	
	M1=allocation(ligne_M1,colonne_M1);
	mat_zeros(M1, ligne_M1,colonne_M1);

	
	//int **vect_alea=runif(0,5,5);
	int *vect=(int *)malloc(1000* sizeof(int));
	for (int i = 0; i < 1000; ++i)
	{
		if (i<500)
		{*vect=-i;}
		else
		{
			*vect=i;
		}
	}
	
	int **vect_alea=rnorm(10,100,vect,1000);
	affichage(vect_alea,1000,1);

//----------------------------------------------------------------


	printf("***Affichage de M1:***\n");
	affichage(M1,ligne_M1,colonne_M1);


 //liberation mémoire
	
	free(vect);
	liberation(vect_alea,1000);
	liberation(M1,ligne_M1);


	return 0;
}