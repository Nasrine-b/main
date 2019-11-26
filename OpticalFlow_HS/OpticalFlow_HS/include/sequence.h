/* ------------------ */
/* --- sequence.h --- */
/* ------------------ */

/*
* Copyright (c) 1999-2009 Lionel Lacassagne, IEF, all rights reserved
 */

#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#ifdef VERBOSE_PRAGMA
//#pragma message ("- *** include sequence.h ***")
#endif

#ifdef __cplusplus
#ifdef PRAGMA_VERBOSE
#pragma message ("C++")
#endif
extern "C" {
#endif

typedef struct {
  char *src_path;
  char *filename;
  char *dst_path;

  int ndigit;          // nombre de digit pour indicer les fichier
  //teFileType filetype; // extension du fichier
  
  int height;         // parametres des images saisie sur la ligne de commande           */
  int width;
  int t;              // time = numero de l'image courante                               */
  int len;            // longueur de la sequence                                         */
  int tstart;         // indice de debut de la sequence                         */
  int tstop;          // indice de fin de la sequence                           */
  int tstep;          // pour sauter des images                                          */
  
  BOOL save;          // save processed files                                            */
  BOOL save_debug;    // save debug file, for visualization                              */
  
} tsSequence;

void        Sequence_Constructor(tsSequence **sequence);
tsSequence* Sequence_pConstructor(void);

void Sequence_Destructor(tsSequence **sequence);
void Sequence_pDestructor(tsSequence *sequence);

void Sequence_Set_SrcPath  (tsSequence *sequence, char *src_path);
void Sequence_Set_Filename (tsSequence *sequence, char *filename);
void Sequence_Set_DstPath  (tsSequence *sequence, char *dst_path);
void Sequence_Set_NDigit   (tsSequence *sequence, int ndigit);
//void Sequence_Set_Filetype (tsSequence *sequence, teFileType filetype);

void Sequence_Set_Dimension(tsSequence *sequence, int hauteur, int largeur);

void Sequence_Set_Width    (tsSequence *sequence, int width);
void Sequence_Set_Height   (tsSequence *sequence, int height);

void Sequence_Set_Largeur    (tsSequence *sequence, int largeur);
void Sequence_Set_Height   (tsSequence *sequence, int hauteur);

void Sequence_Set_Len      (tsSequence *sequence, int len);

void Sequence_Set_T        (tsSequence *sequence, int t);
void Sequence_Set_TStart   (tsSequence *sequence, int tstart);
void Sequence_Set_TStop    (tsSequence *sequence, int tstop);
void Sequence_Set_TStep    (tsSequence *sequence, int tstep);
void Sequence_Set_Len      (tsSequence *sequence, int len);

void Sequence_Set_Save     (tsSequence *sequence, BOOL save);
void Sequence_Set_SaveDebug(tsSequence *sequence, BOOL save_debug);

void Sequence_Parse_File           (char *filename, tsSequence  *sequence);
void Sequence_Display_Configuration(tsSequence *sequence);

void Sequence_Save_Configuration   (tsSequence *sequence, char *filename_cfg);
void Sequence_Load_Configuration   (char *filename_cfg, tsSequence *sequence);

void Sequence_Next(tsSequence *sequence);

char*Sequence_Get_CompleteFilename(tsSequence *sequence);

char* Sequence_Get_SrcPath (tsSequence *sequence);
char* Sequence_Get_Filename(tsSequence *sequence);
char* Sequence_Get_DstPath (tsSequence *sequence);
int   Sequence_Get_NDigit  (tsSequence *sequence);

int Sequence_Get_T     (tsSequence *sequence);
int Sequence_Get_TStart(tsSequence *sequence);
int Sequence_Get_TStop (tsSequence *sequence);
int Sequence_Get_TStep (tsSequence *sequence);
int Sequence_Get_Len   (tsSequence *sequence);

int Sequence_Get_Height(tsSequence *sequence);
int Sequence_Get_Width (tsSequence *sequence);

BOOL Sequence_Get_Save     (tsSequence *sequence);
BOOL Sequence_Get_SaveDebug(tsSequence *sequence);

#ifdef __cplusplus
}
#endif

#endif // __SEQUENCE_H__
