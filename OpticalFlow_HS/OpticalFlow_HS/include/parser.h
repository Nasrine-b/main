/* ---------------- */
/* --- parser.c --- */
/* ---------------- */

/*
 * Basic Command Line Parser
 */

/*
 * Copyright (c) 1999-2009 Lionel Lacassagne
 *
 * 1999-04-12
 * 2002-02-13 : ajout du type BOOL
 */

#ifndef __PARSER_H__
#define __PARSER_H__

#ifdef VERBOSE_PRAGMA
//#pragma message ("- *** include parser.h ***")
#endif

#ifdef __cplusplus
#ifdef PRAGMA_VERBOSE
#pragma message ("C++")
#endif
extern "C" {
#endif

#define scanLINE "%81[^\n]%*[\n]"
#define scanINT  "*[ \t]%d"
#define scanCHAR "*[ \t]%1s"
#define scanSTR0 "%*[^ ^\t]%81s"
#define scanSTR  "*81[^ ^\t]%81s"

void Parser_ScanLine (FILE *file, char *line);
void Parser_ScanLine2(FILE *file, char *line);
void Parser_ScanChar (char *line, char *word, char  *c);
void Parser_ScanStr  (char *line, char *word, char  *str);
void Parser_ScanInt  (char *line, char *word, int   *integer);
void Parser_ScanFloat(char *line, char *word, float *real);
void Parser_ScanBOOL (char *line, char *word, BOOL  *b);

void Parser_Reset     (void);
void Parser_WriteChar (char *word, char  c);
void Parser_WriteStr  (char *word, char *str);
void Parser_WriteInt  (char *word, int   integer);
void Parser_WriteFloat(char *word, float real);
void Parser_WriteBOOL (char *word, BOOL  b);

#ifdef __cplusplus
}
#endif

#endif // __PARSER_H__
