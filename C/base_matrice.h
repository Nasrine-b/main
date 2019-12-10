#ifndef __CREATION_H_
#define __CREATION_H_

int ** allocation(int l, int c);
void liberation(int **M, int l);
void mat_zeros(int **M, int nb_ligne, int nb_colonne);
void affichage(int **M, int nb_ligne, int nb_colonne);

#endif