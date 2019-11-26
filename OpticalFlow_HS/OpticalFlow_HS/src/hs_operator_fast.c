/* ----------------------------- */
/* --- HornShunck_operator.c --- */
/* ----------------------------- */

/*
 * Copyright (c) 2017, Lionel Lacassagne, all rights reserved
 * LIP6, UPMC, CNRS
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// NRC lib
#include "nrc.h"
#include "nrdef.h"
#include "nrtype.h"
#include "nralloc.h"
#include "nrarith.h"
#include "nrio.h"

// local
//#include "parser.h"
//#include "sequence.h"

//#include "of_config.h"
//#include "of.h"
#include "of_macro.h"
#include "of_function.h"

//#include "hs_algo.h"
#include "hs_operator.h"

#include "interp.h"
#include "interp_extra.h"

#ifdef OPENMP
#include <omp.h>
#endif

