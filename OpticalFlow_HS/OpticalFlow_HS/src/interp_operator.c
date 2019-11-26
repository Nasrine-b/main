/* ------------------------- */
/* --- interp_operator.c --- */
/* ------------------------- */

/*
* interpolation functions
*/

/*
 * Copyright (c) 2017-2017 Lionel Lacassagne, all rights reserved, LIP6, UPMC, CNRS
 * Copyright (c) 2018-2018 Lionel Lacassagne, all rights reserved, LIP6, Sorbonne Universite, CNRS
 */

/*
 * 2017: creation des fonctions d'interpolation cubique
 * 2018: ajout de la fonction d'interpolation lineaire
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "nrc.h"
#include "nrdef.h"
#include "nrtype.h"

#include "of_macro.h"
#include "interp_operator.h"

// ------------------------------------------------------------
float32 linear_interpolation(float32 f1, float32 f2, float32 x)
// ------------------------------------------------------------
{
    float32 c1 = 1.0f - x;
    float32 c2 = x;
    
    float32 f = c1 * f1 + c2 * f2;
    
    return f;
}
// ----------------------------------------------------------------------------------------
float32 cubic_interpolation_poly(float32 f1, float32 f2, float32 f3, float32 f4, float32 x)
// ----------------------------------------------------------------------------------------
{
    float32 x2 = x * x;
    float32 x3 = x * x2;
    
    float32 c1 = -0.5f * (    x3 - 2 * x2 + x);
    float32 c2 =  0.5f * (3 * x3 - 5 * x2 + 2);
    float32 c3 = -0.5f * (3 * x3 - 4 * x2 - x);
    float32 c4 =  0.5f * (    x3 -     x2    );
    
    float32 f = c1 * f1 + c2 * f2 + c3 * f3 + c4 * f4;
    
    return f;
}
// ------------------------------------------------------------------------------------------
float32 cubic_interpolation_horner(float32 f1, float32 f2, float32 f3, float32 f4, float32 x)
// ------------------------------------------------------------------------------------------
{
    float32 x2 = x * x;
    
    float32 c1 = - ((x - 2) * x + 1) * x;
    float32 c2 = (3 * x - 5) * x2 + 2;
    float32 c3 = - ((3 * x - 4) * x  - 1) * x;
    float32 c4 = (x - 1) * x2;
    
    float32 f = 0.5f * (c1 * f1 + c2 * f2 + c3 * f3 + c4 * f4);
    
    return f;
}
// ----------------------------------------------------------------------------------------
float32 cubic_interpolation_cell(float32 f1, float32 f2, float32 f3, float32 f4, float32 x)
// ----------------------------------------------------------------------------------------
{
    // from IPOL
    
    float32 f = f2 + 0.5 * x * (f3 - f1 + x * (2 *  f1 - 5 * f2 + 4 * f3 - f4 + x * (3 * (f2 - f3) + f4 - f1)));
    return f;
}
// ----------------------------------------------------------------------------------------
float32 cubic_interpolation_ipol(float32 f1, float32 f2, float32 f3, float32 f4, float32 x)
// ----------------------------------------------------------------------------------------
{
    // from IPOL
    float32 c3 = 3 * (f2 - f3) + f4 - f1;
    float32 c2 = 2 *  f1 - 5 * f2 + 4 * f3 - f4;
    float32 c1 = f3 - f1;
    float32 c0 = f2;
    
    /*printf("c3 = %6.2f\n", c3);
    printf("c2 = %6.2f\n", c2);
    printf("c1 = %6.2f\n", c1);
    printf("c0 = %6.2f\n", c0);*/
    printf("%8.2f%8.2f%8.2f%8.2f\n", c3, c2, c1, c0);
    
    //float32 f = f2 + 0.5 * x * (f3 - f1 + x * (2 *  f1 - 5 * f2 + 4 * f3 - f4 + x * (3 * (f2 - f3) + f4 - f1)));
    float32 f = ((c3 * x + c2) * x + c1) * 0.5f * x + c0;
    return f;
}
// -----------------------------------------------------------------------------------------
float32 cubic_interpolation_burke(float32 f0, float32 f1, float32 f2, float32 f3, float32 x)
// -----------------------------------------------------------------------------------------
{
    // from Paul Burke
    // (http://paulbourke.net/miscellaneous/interpolation/)

    // formule "fausse" car interpolation exacte que aux points f1 (x=0) et f2 (x=1)
    
    float32 x2 = x * x;
    float32 x3 = x * x2; // only for poly
    
    float32 c0 = f3 - f2 - f0 + f1;
    float32 c1 = f0 - f1 - c0;
    float32 c2 = f2 - f0;
    float32 c3 = f1;
    
    return(c0 * x3 + c1 *x2 + c2 * x + c3);
    
    //float32 fp = c1 * x3 + c2 * x2 + c3 * x + c4;   // poly
    //float32 fh = ((c1 * x + c2) * x + c3) * x + c4; // horner
    //return fp;
}
// ---------------------------------------------------------------------------------------------
float32 cubic_interpolation_breeuwsma(float32 f1, float32 f2, float32 f3, float32 f4, float32 x)
// ---------------------------------------------------------------------------------------------
{
    // from Paul Breeuwsma also known as Catmull-Rom splines
    // (http://paulbourke.net/miscellaneous/interpolation/)
    
    float32 x2 = x * x;
    float32 x3 = x * x2;
    
    float32 c1 = -0.5 * f1 + 1.5 * f2 - 1.5 * f3 + 0.5 * f4;
    float32 c2 = f1 - 2.5 * f2 + 2 * f3 - 0.5 * f4;
    float32 c3 = -0.5 * f1 + 0.5 * f3;
    float32 c4 = f2;
    
    float32 fp = c1 * x3 + c2 * x2 + c3 * x + c4;   // poly
    //float32 fh = ((c1 * x + c2) * x + c3) * x + c4; // horner
    
    return fp;
}
