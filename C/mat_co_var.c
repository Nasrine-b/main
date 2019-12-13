#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "base_matrice.h"
#include "mat_co_var.h"


int ****allocation_sous_mat(int dim1, int dim2)
{
	
	int ****sous_M;
	sous_M = (int ****) malloc( dim1* sizeof(int ***) );
	for(int i=0; i < dim1; i++)
	{
		sous_M[i]=(int ***)malloc(dim2* sizeof(int **));
	
	}

	for (int i = 0; i < dim1; ++i)
	{
		for (int j = 0; j < dim2; ++j)
		{
			sous_M[i][j] = allocation(2, 2);
		}
		
	}
	return sous_M;

}
void sous_mat_zeros(int ****sous_M, int dim1, int dim2)
{
	
	for (int i = 0; i < dim1; ++i)
	{
		for (int j = 0; j < dim2; ++j)
		{
			mat_zeros(sous_M[i][j], 2, 2);
		}
		
	}

}

void liberation_sous_mat(int ****res, int dim1, int dim2)
{
	
	for (int i = 0; i < dim1; ++i)
	{
		for (int j = 0; j < dim2; ++j)
		{
			liberation(res[i][j], 2);
		}		
	}
	for (int i = 0; i < dim1; ++i)
	{
		free(res[i]);
	}
	free(res);

}
void affichage_sous_mat(int ****res, int dim1, int dim2)
{
	for (int i = 0; i < dim1; ++i)
	{
		for (int j = 0; j < dim2; ++j)
		{
			printf("Mat[%d][%d]\n",i,j );
			affichage(res[i][j],2,2);
		}
	}
}

// int ****mat_cov2(int **vect, int dim1,int dim2)
// {

// 	int ****res=allocation_sous_mat(dim1,dim2);
// 	sous_mat_zeros(res,dim1, dim2);


//  	for (unsigned int i=0 ; i < dim1; i++)
//  	{
//  	    for (unsigned int j=0; j < dim2; j++)
//  	    {


//  	        res(:,:,i,j)=cov_2_2_mat(Vec(i,:),Vec(j,:));
//  	    }
 
//  	}
// } 

//function [matrix] = cov_2_2_mat(X1,X2)
int **cov_2_2_mat(int *vec1,int *vec2)
{
	int **res=allocation(2,2);
	mat_zeros(res,2,2);

	int elem1=somme_elem_mat(vec1,2)/2;
	int elem2=somme_elem_mat(vec2,2)/2;

	res[0][0]=(vec1[0]-elem1)*(vec1[0]-elem1)+(vec1[1]-elem1)*(vec1[1]-elem1);
	res[0][0]=(vec2[0]-elem2)*(vec2[0]-elem2)+(vec2[1]-elem2)*(vec2[1]-elem2);
	res[0][1]=(vec1[0]-elem1)*(vec2[0]-elem2)+(vec1[1]-elem1)*(vec2[1]-elem2);
	res[1][0]=res[0][1];
	return res;
}
int somme_elem_mat(int *M, int l)
{
	int res;
	for (int i = 0; i < l; ++i)
	{
		
		res=M[i]+res;
		
	}
	return res;

 }

// matrix=zeros(2,2);
// x1_m=sum(X1)/2;
// x2_m=sum(X2)/2;
// matrix(1,1)=(X1(1)-x1_m)^2+(X1(2)-x1_m)^2;
// matrix(2,2)=(X2(1)-x2_m)^2+(X2(2)-x2_m)^2;
// matrix(1,2)=(X1(1)-x1_m)*(X2(1)-x2_m)+(X1(2)-x1_m)*(X2(2)-x2_m);
// matrix(2,1)=matrix(1,2);

// end