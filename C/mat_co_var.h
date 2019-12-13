#ifndef __mat_co_var_H_
#define __mat_co_var_H_


int ****allocation_sous_mat(int dim1, int dim2);
void sous_mat_zeros(int ****sous_M, int dim1, int dim2);
void liberation_sous_mat(int ****res, int dim1, int dim2);
void affichage_sous_mat(int ****res, int dim1, int dim2);
int **cov_2_2_mat(int *vec1,int *vec2);
int somme_elem_mat(int *M, int l);
// int ****mat_cov2(int **vect, int dim1,int dim2);
#endif