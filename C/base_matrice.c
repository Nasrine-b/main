#include <stdio.h>
#include <stdlib.h>

int ** allocation(int l, int c)
{
	int **M;
	//Allocation mémoire
	M = (int **) malloc( l* sizeof(int *) );
	//vérification allocation mémoire
	if (M == NULL)
	{
		printf("Erreur allocation mémoire\n");
		exit(1);
	}
	for(int i=0; i<l;i++)
	{
		M[i]=(int *)malloc(c* sizeof(int));
		//vérification allocation mémoire
		if (M[i] == NULL)
		{
			printf("Erreur allocation mémoire\n");
			exit(1);
		}
	}
	return M;
}

void liberation(int **M, int l)
{
	//liberation mémoire
	for(int i=0; i<l;i++)
	{
		free(M[i]);
	}
	free(M);
}
/*fonction innitialise une matrice
*param M une matrice pointeur
*param nb_ligne entier 
*param nb_colonne un entier
*return initialise la matrice avec des coefs choisis par l'utilisateur 
*/
void mat_zeros(int **M, int nb_ligne, int nb_colonne)
{
	for (int i = 0; i < nb_ligne; i++)
	{
		for (int j = 0; j < nb_colonne; j++)
		{
			M[i][j]=0;
		}
	}
}
/*fonction affiche une matrice
*param M une matrice pointeur
*param nb_ligne entier 
*param nb_colonne un entier
*return affichage 
*/
void affichage(int **M, int nb_ligne, int nb_colonne)
{
	for (int i = 0; i < nb_ligne; i++)
	{
		for (int j = 0; j < nb_colonne; j++)
		{
			printf(" %d ", M[i][j]);
		}
		printf("\n");
	}
}