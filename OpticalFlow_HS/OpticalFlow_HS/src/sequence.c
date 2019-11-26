/* ------------------- */
/* --- sequence.c ---- */
/* ------------------- */

/*
* Copyright (c) 2006-2009 Lionel Lacassagne, IEF, all rights reserved
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h> // va_list
#include <string.h>

// NRC lib
#include "nrc.h"
#include "nrdef.h"
#include "nrtype.h"
#include "nralloc.h"
#include "nrarith.h"
#include "nrio.h"

#include "parser.h"

//#include "img.h" // pour champs ximage
#include "sequence.h"

// -----------------------------------
void Sequence_Error(char *format, ...)
// -----------------------------------
/* global error handler */
{
    va_list args;
    FILE *file;
    
    va_start(args, format);
    
    file = fopen("error.txt", "wb");
    if(file == NULL) {
        fprintf (stderr,"Sequence_Error : can't write error to disk\n");
    } else {
        fprintf (file,"Sequence_Error\n");
        vfprintf(file, format, args);
        fprintf (file, "\n");
        fclose  (file);
    }
    fprintf (stderr, "Sequence_Erreur : ");
    vfprintf(stderr, format, args);
    fprintf (stderr, "\n");
    
    va_end(args);
    
    /*getchar();*/

  exit(-1);
}
// ---------------------------------------------
void Sequence_Constructor(tsSequence **sequence)
// ---------------------------------------------
{
    *sequence = Sequence_pConstructor();
}
// ------------------------------------
tsSequence* Sequence_pConstructor(void)
// ------------------------------------
{
    tsSequence *sequence;
    
    sequence = (tsSequence*) malloc(sizeof(tsSequence));
    if(sequence == NULL) {
        Sequence_Error("Sequence_pConstructor error\n");
        return NULL;
    }
    sequence->src_path = (char*) malloc(1024*sizeof(char));
    sequence->filename = (char*) malloc(1024*sizeof(char));
    sequence->dst_path = (char*) malloc(1024*sizeof(char));
    
    if( (sequence->src_path==NULL) || (sequence->filename==NULL) || (sequence->dst_path==NULL))
        Sequence_Error("Sequence_pConstructor malloc error\n");
    
    return sequence;
}
// --------------------------------------------
void Sequence_Destructor(tsSequence **sequence)
// --------------------------------------------
{
    Sequence_pDestructor(*sequence);
}
// --------------------------------------------
void Sequence_pDestructor(tsSequence *sequence)
// --------------------------------------------
{
    free(sequence->src_path); sequence->src_path = NULL;
    free(sequence->filename); sequence->filename = NULL;
    free(sequence->dst_path); sequence->dst_path = NULL;
    
    free(sequence);
}
// ------------------------------------------------------------
void Sequence_Set_SrcPath(tsSequence *sequence, char *src_path)
// ------------------------------------------------------------
{
    strcpy(sequence->src_path, src_path);
}
// -------------------------------------------------------------
void Sequence_Set_Filename(tsSequence *sequence, char *filename)
// -------------------------------------------------------------
{
    strcpy(sequence->filename, filename);
}
// ------------------------------------------------------------
void Sequence_Set_DstPath(tsSequence *sequence, char *dst_path)
// ------------------------------------------------------------
{
    strcpy(sequence->dst_path, dst_path);
}
// -------------------------------------------------------
void Sequence_Set_NDigit(tsSequence *sequence, int ndigit)
// -------------------------------------------------------
{
    sequence->ndigit = ndigit;
}
// --------------------------------------------------------------------
//void Sequence_Set_Filetype(tsSequence *sequence, teFileType filetype)
// --------------------------------------------------------------------
/*{
    sequence->filetype = filetype;
}*/
// ---------------------------------------------------------------------
void Sequence_Set_Dimension(tsSequence *sequence, int height, int width)
// ---------------------------------------------------------------------
{
    sequence->height = height;
    sequence->width = width;
}
// -----------------------------------------------------
void Sequence_Set_Width(tsSequence *sequence, int width)
// -----------------------------------------------------
{
    sequence->width = width;
}
// -------------------------------------------------------
void Sequence_Set_Height(tsSequence *sequence, int height)
// -------------------------------------------------------
{
    sequence->height = height;
}
// -------------------------------------------------------
void Sequence_Set_Largeur(tsSequence *sequence, int width)
// -------------------------------------------------------
{
    sequence->width = width;
}
// --------------------------------------------------------
void Sequence_Set_Hauteur(tsSequence *sequence, int height)
// --------------------------------------------------------
{
    sequence->height = height;
}
// -------------------------------------------------
void Sequence_Set_Len(tsSequence *sequence, int len)
// -------------------------------------------------
{
    sequence->len = len;
}
// ---------------------------------------------
void Sequence_Set_T(tsSequence *sequence, int t)
// ---------------------------------------------
{
    sequence->t = t;
}
// -------------------------------------------------------
void Sequence_Set_TStart(tsSequence *sequence, int tstart)
// -------------------------------------------------------
{
    sequence->tstart = tstart;
    sequence->t      = tstart;
}
// -----------------------------------------------------
void Sequence_Set_TStop(tsSequence *sequence, int tstop)
// -----------------------------------------------------
{
    sequence->tstop = tstop;
}
// -----------------------------------------------------
void Sequence_Set_TStep(tsSequence *sequence, int tstep)
// -----------------------------------------------------
{
    sequence->tstep = tstep;
}
// ----------------------------------------------------
void Sequence_Set_Save(tsSequence *sequence, BOOL save)
// ----------------------------------------------------
{
    sequence->save = save;
}
// ---------------------------------------------------------------
void Sequence_Set_SaveDebug(tsSequence *sequence, BOOL save_debug)
// ---------------------------------------------------------------
{
    sequence->save_debug = save_debug;
}
// ---------------------------------------------------------------
void Sequence_Parse_File(char *filename_cfg, tsSequence *sequence)
// ---------------------------------------------------------------
{
    char src_path[1024];
    char filename[80];
    char dst_path[1024];
    
    char line[1024];
    //char str[1024];
    
    int height, width;
    int len;
    int tstart, tstop, tstep;
    int ndigit;
    //teFileType filetype;
    
    BOOL save, save_debug;
    
    FILE *file;
    
    // --- First step: read parameters ----
    
    file = fopen(filename_cfg, "rb");
    if(file == NULL) {
        Sequence_Error("Sequence_Parse_File : can't open %s\n", filename_cfg);
    }
    
    Parser_Reset();
    
    Parser_ScanLine2(file, line); Parser_ScanStr(line, "src_path", src_path);
    Parser_ScanLine2(file, line); Parser_ScanStr(line, "filename", filename);
    Parser_ScanLine2(file, line); Parser_ScanInt(line, "ndigit", &ndigit);
    //Parser_ScanLine2(file, line); Parser_ScanStr(line, "filetype", str);
    
    /*filetype = NO_FILETYPE;
    if( strcmp(str, pcFileType[FILETYPE_PGM]) == 0 ) filetype = FILETYPE_PGM;
    if( strcmp(str, pcFileType[FILETYPE_PNG]) == 0 ) filetype = FILETYPE_PNG;
    if(!filetype) Sequence_Error("Bad filetype : %s", str);*/
    
    Parser_ScanLine2(file, line); Parser_ScanStr(line, "dst_path", dst_path);
    Parser_ScanLine2(file, line); Parser_ScanInt(line, "width", &width);
    Parser_ScanLine2(file, line); Parser_ScanInt(line, "height", &height);
    Parser_ScanLine2(file, line); Parser_ScanInt(line, "tstart", &tstart);
    Parser_ScanLine2(file, line); Parser_ScanInt(line, "tstop", &tstop);
    Parser_ScanLine2(file, line); Parser_ScanInt(line, "tstep", &tstep);
    
    Parser_ScanLine2(file, line); Parser_ScanBOOL(line, "save", &save);
    Parser_ScanLine2(file, line); Parser_ScanBOOL(line, "save_debug", &save_debug);
    
    len = (tstop - tstart + 1) / tstep;

    Sequence_Set_SrcPath  (sequence, src_path);
    Sequence_Set_Filename (sequence, filename);
    Sequence_Set_NDigit   (sequence, ndigit);
    //Sequence_Set_Filetype (sequence, filetype);
    Sequence_Set_DstPath  (sequence, dst_path);
    
    Sequence_Set_Width    (sequence, width);
    Sequence_Set_Height   (sequence, height);
    
    Sequence_Set_T        (sequence, tstart);
    Sequence_Set_TStart   (sequence, tstart);
    Sequence_Set_TStop    (sequence, tstop);
    Sequence_Set_TStep    (sequence, tstep);
    Sequence_Set_Len      (sequence, len);

    Sequence_Set_Save     (sequence, save);
    Sequence_Set_SaveDebug(sequence, save_debug);
    
    }
// ------------------------------------------------------
void Sequence_Display_Configuration(tsSequence *sequence)
// ------------------------------------------------------
{
    char *src_path = sequence->src_path;
    char *filename = sequence->filename;
    char *dst_path = sequence->dst_path;
    
    int height = sequence->height;
    int width = sequence->width;
    //int len     = sequence->len;
    int tstart  = sequence->tstart;
    int tstop   = sequence->tstop;
    //int tstep   = sequence->tstep;
    int ndigit  = sequence->ndigit;
    
    //teFileType filetype = sequence->filetype;
    
    //printf("%s %s.%s -> %s\n", src_path, filename, pcFileType[sequence->filetype], dst_path);
    printf("%s %s -> %s\n", src_path, filename, dst_path);
    printf("%d x %d [%d..%d].%d\n", width, height, tstart, tstop, ndigit);
}
// -----------------------------------------------------------------------
void Sequence_Save_Configuration(tsSequence *sequence, char *filename_cfg)
// -----------------------------------------------------------------------
{
    FILE *f;
    
    f = fopen(filename_cfg, "rb");
    if(f==NULL) {
        Sequence_Error("Sequence_Save_Configuration: can' open file %s\n", filename_cfg);
    }
    
    fprintf(f, "src_path = %s\n", sequence->src_path);
    fprintf(f, "filename = %s\n", sequence->filename);
    fprintf(f, "dst_path = %s\n", sequence->dst_path);
    fprintf(f, "width  = %d\n", sequence->width);
    fprintf(f, "height  = %d\n", sequence->height);
    fprintf(f, "tstart   = %d\n", sequence->tstart);
    fprintf(f, "tstop    = %d\n", sequence->tstop);
    fprintf(f, "tstep    = %d\n", sequence->tstep);
    fprintf(f, "ndigit   = %d\n", sequence->ndigit);
    //fprintf(f, "filetype = %s\n", pcFileType[sequence->filetype]);
    
}
// -----------------------------------------------------------------------
void Sequence_Load_Configuration(char *filename_cfg, tsSequence *sequence)
// -----------------------------------------------------------------------
{
    Sequence_Parse_File(filename_cfg, sequence);
}
// -------------------------------------
void Sequence_Next(tsSequence *sequence)
// -------------------------------------
{
    sequence->t += sequence->tstep;
}
// ------------------------------------------------------
char* Sequence_Get_CompleteFilename(tsSequence *sequence)
// ------------------------------------------------------
{
    static char complete_filename[1024];
    
    char *src_path = sequence->src_path;
    char *filename = sequence->filename;
    int t = sequence->t;
    int ndigit = sequence->ndigit;
    
    generate_path_filename_k_ndigit_extension  (src_path, filename, t, ndigit, "pgm", complete_filename);
    
    return complete_filename;
}
// -------------------------------------
int Sequence_Get_T(tsSequence *sequence)
// -------------------------------------
{
    return sequence->t;
}
// ------------------------------------------
int Sequence_Get_TStart(tsSequence *sequence)
// ------------------------------------------
{
    return sequence->tstart;
}
// -----------------------------------------
int Sequence_Get_TStop(tsSequence *sequence)
// -----------------------------------------
{
    return sequence->tstop;
}
// -----------------------------------------
int Sequence_Get_TStep(tsSequence *sequence)
// -----------------------------------------
{
    return sequence->tstep;
}
// ---------------------------------------
int Sequence_Get_Len(tsSequence *sequence)
// ---------------------------------------
{
    return sequence->len;
}
// -----------------------------------------
int Sequence_Get_Width(tsSequence *sequence)
// -----------------------------------------
{
    return sequence->width;
}
// ------------------------------------------
int Sequence_Get_Height(tsSequence *sequence)
// ------------------------------------------
{
    return sequence->height;
}
// ---------------------------------------------
char* Sequence_Get_SrcPath(tsSequence *sequence)
// ---------------------------------------------
{
    return sequence->src_path;
}
// ----------------------------------------------
char* Sequence_Get_Filename(tsSequence *sequence)
// ----------------------------------------------
{
    return sequence->filename;
}
// ---------------------------------------------
char* Sequence_Get_DstPath(tsSequence *sequence)
// ---------------------------------------------
{
    return sequence->dst_path;
}
// ------------------------------------------
int Sequence_Get_NDigit(tsSequence *sequence)
// ------------------------------------------
{
    return sequence->ndigit;
}
// -----------------------------------------
BOOL Sequence_Get_Save(tsSequence *sequence)
// -----------------------------------------
{
    return sequence->save;
}
// ----------------------------------------------
BOOL Sequence_Get_SaveDebug(tsSequence *sequence)
// ----------------------------------------------
{
    return sequence->save_debug;
}
