/* -------------------------- */
/* --- test_hs_operator.c --- */
/* -------------------------- */

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

#include "nrpyramid.h"

// local
#include "of_macro.h"

#include "of_function.h"
//#include "of_function_fix.h"

#include "hs_operator.h"
//#include "hs_operator_fix.h"
//#include "hs_operator_fast.h"

#include "interp.h"
//#include "interp_fix.h"

#include "interp_operator.h"
//#include "interp_operator_fix.h"

#include "of_score_op.h"

// ---------------------------------------------------------------------------------------------------------
void display_si64matrix_quantif(sint64 **X, int i0, int i1, int j0, int j1, int q, char *format, char *name)
// ---------------------------------------------------------------------------------------------------------
{
    if(name) puts(name);

    float32 Q = (float32) CALC_Q(q);

    for(int i= i0; i <= i1; i++) {
        for(int j = j0; j <= j1; j++) {
            float32 x = (float32) X[i][j];
            float32 f = x / Q;
            printf(format, f);
        }
        putchar('\n');
    }
    putchar('\n');

}
// ------------------------------------------
int HS_test_mono1_F32(int argc, char *argv[])
// ------------------------------------------
{
    int height, width;  // image size
    int i0, i1, j0, j1; // image intervals
    int b = 4;          // border

    float32 alpha2;

    double ae, ee;                 // angle and endpoint errors
    double ie_nir, ie_lin, ie_cub; // interpolation error for {nearest, linear, cubic} interpolation

    char *src_path = "/nfs/home/sasl/eleves/main/3800549/Bureau/OpticalFlow_HS/OpticalFlow_HS/config/";
    char *dst_path = "/nfs/home/sasl/eleves/main/3800549/Bureau/OpticalFlow_HS/OpticalFlow_HS/config/";

    char *filename0;
    char *filename1;

    char complete_filename [1024];
    char complete_filename0[1024];
    char complete_filename1[1024];

    uint8 **I0;
    uint8 **I1;

    float32 **E0, **E1;       // pair of images (t), (t+1)
    float32 **Ex, **Ey, **Et; // {x, y, t} gradient

    float32 **U0, **V0; // input flow
    float32 **U1, **V1; // output flow
    float32 **Um, **Vm; // average flow (moyen)

    float32 **Ugt, **Vgt; // ground-true

    float32 **Er_nir, **Er_lin, **Er_cub;
    uint8   **Ir_nir, **Ir_lin, **Ir_cub;

    puts("-------------------------");
    puts("--- HS_test_mono1_F32 ---");
    puts("-------------------------");

    filename0 = "I0.pgm"; filename1 = "I1.pgm"; //synthetic cosine/sine
    //filename0 = "coastguard060.pgm"; filename1 = "coastguard061.pgm";

    generate_path_filename(src_path, filename0, complete_filename0);
    generate_path_filename(src_path, filename1, complete_filename1);

    // load & alloc
    I0 = LoadPGM_ui8matrix(complete_filename0, &i0, &i1, &j0, &j1);
    I1 = LoadPGM_ui8matrix(complete_filename1, &i0, &i1, &j0, &j1);

    height = i1 - i0 + 1;
    width  = j1 - j0 + 1;
    printf("size = %d x %d\n", width, height);

    // alloc
    E0 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    E1 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Ex = f32matrix(i0, i1, j0, j1);
    Ey = f32matrix(i0, i1, j0, j1);
    Et = f32matrix(i0, i1, j0, j1);

    U0 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    V0 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    U1 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    V1 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Um = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    Vm = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Ugt = f32matrix(i0, i1, j0, j1);
    Vgt = f32matrix(i0, i1, j0, j1);

    Er_nir = f32matrix(i0, i1, j0, j1);
    Er_lin = f32matrix(i0, i1, j0, j1);
    Er_cub = f32matrix(i0, i1, j0, j1);

    Ir_nir = ui8matrix(i0, i1, j0, j1);
    Ir_lin = ui8matrix(i0, i1, j0, j1);
    Ir_cub = ui8matrix(i0, i1, j0, j1);

    // init
    zero_f32matrix(Ex, i0, i1, j0, j1);
    zero_f32matrix(Ey, i0, i1, j0, j1);
    zero_f32matrix(Et, i0, i1, j0, j1);

    zero_f32matrix(E0, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_f32matrix(E1, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_f32matrix(U0, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_f32matrix(V0, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_f32matrix(U1, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_f32matrix(V1, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_f32matrix(Um, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_f32matrix(Vm, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_f32matrix(Ugt, i0, i1, j0, j1);
    zero_f32matrix(Vgt, i0, i1, j0, j1);

    set_f32matrix(Ugt, i0, i1, j0, j1, 0.75f);
    set_f32matrix(Vgt, i0, i1, j0, j1, 0.25f);

    // convert = division par 255
    convert_ui8matrix_f32matrix_quantif(I0, i0, i1, j0, j1, 255, E0);  // attention: pas 256
    convert_ui8matrix_f32matrix_quantif(I1, i0, i1, j0, j1, 255, E1);

    //display_f32matrix(E0, height/4, 3*height/4, width/4, 3*width/4, "%8.2f", "E0");
    //display_f32matrix(E1, height/4, 3*height/4, width/4, 3*width/4, "%8.2f", "E1");

    // calc

    alpha2 = 0.635f; // mauvaise convergence
    alpha2 = 0.05f;  // bonne convergence
    //alpha2 = 0.04f;

    hs_step0_F32(E0, E1, height, width, alpha2, Ex, Ey, Et, U1, V1);

    //display_f32matrix(Ex, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "Ex");
    //display_f32matrix(Ey, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "Ey");
    //display_f32matrix(Et, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "Et");

    ae = 0.0; ee = 0.0;
    for(int iter=0; iter<500; iter++) {

        hs_1step_F32(U0, V0, Ex, Ey, Et, height, width, alpha2, Um, Vm, U1, V1);
        ae = of_calcAngularError_F32 (U0, V0, Ugt, Vgt, height, width);
        ee = of_calcEndpointError_F32(U0, V0, Ugt, Vgt, height, width);
        printf("iter %3d angle = %8.3f  norm = %8.3f \n", iter, ae, ee);
        SWAP_F32(U0, U1);
        SWAP_F32(V0, V1);
    }
    display_f32matrix(U0, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "U");
    display_f32matrix(V0, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "V");

    generate_path_filename(dst_path, "U_.txt", complete_filename);  write_f32matrix(U0, i0, i1, j0, j1, "%8.3f", complete_filename);
    generate_path_filename(dst_path, "V_.txt", complete_filename);  write_f32matrix(V0, i0, i1, j0, j1, "%8.3f", complete_filename);

    interpNearest_F32(E1, U0, V0, height, width, Er_nir);
    interpLinear_F32 (E1, U0, V0, height, width, Er_lin);
    interpCubic_F32  (E1, U0, V0, height, width, Er_cub);

    convert_f32matrix_ui8matrix_quantif(Er_nir, i0, i1, j0, j1, 255, Ir_nir); // attention: pas 256
    convert_f32matrix_ui8matrix_quantif(Er_lin, i0, i1, j0, j1, 255, Ir_lin);
    convert_f32matrix_ui8matrix_quantif(Er_cub, i0, i1, j0, j1, 255, Ir_cub);

    display_f32matrix(Er_cub, height/4, 3*height/4, width/4, 3*width/4, "%6.3f ", "Er_cub");
    display_ui8matrix(Ir_cub, height/4, 3*height/4, width/4, 3*width/4, "%6d ",   "Ir_cub");

    generate_path_filename(dst_path, "Irec_nir_.pgm", complete_filename); SavePGM_ui8matrix(Ir_nir, i0, i1, j0, j1, complete_filename);
    generate_path_filename(dst_path, "Irec_lin_.pgm", complete_filename); SavePGM_ui8matrix(Ir_lin, i0, i1, j0, j1, complete_filename);
    generate_path_filename(dst_path, "Irec_cub_.pgm", complete_filename); SavePGM_ui8matrix(Ir_cub, i0, i1, j0, j1, complete_filename);

    generate_path_filename(dst_path, "Irec_nir_.txt", complete_filename); write_ui8matrix_number(Ir_nir, i0, i1, j0, j1, "%4d", complete_filename);
    generate_path_filename(dst_path, "Irec_lin_.txt", complete_filename); write_ui8matrix_number(Ir_lin, i0, i1, j0, j1, "%4d", complete_filename);
    generate_path_filename(dst_path, "Irec_cub_.txt", complete_filename); write_ui8matrix_number(Ir_cub, i0, i1, j0, j1, "%4d", complete_filename);

    ae = of_calcAngularError_F32 (U0, V0, Ugt, Vgt, height, width);
    ee = of_calcEndpointError_F32(U0, V0, Ugt, Vgt, height, width);

    ie_nir = of_calcInterpolationError_F32(E0, Er_nir, height, width);
    ie_lin = of_calcInterpolationError_F32(E0, Er_lin, height, width);
    ie_cub = of_calcInterpolationError_F32(E0, Er_cub, height, width);

    printf("error angle(deg) = %.3f  error L2norm(pixel) = %.3f\n", ae, ee);
    printf("interp_nir(gray) = %.3f  interp_lin(gray) = %.3f  interp_cub(gray) = %.3f\n", 256*ie_nir, 256*ie_lin, 256*ie_cub);

    // free
    puts("free");
    free_ui8matrix(I0, i0, i1, j0, j1);
    free_ui8matrix(I1, i0, i1, j0, j1);

    free_f32matrix(E0, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(E1, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Ex, i0, i1, j0, j1);
    free_f32matrix(Ey, i0, i1, j0, j1);
    free_f32matrix(Et, i0, i1, j0, j1);

    free_f32matrix(U0, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(V0, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(U1, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(V1, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Um, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(Vm, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Ugt, i0, i1, j0, j1);
    free_f32matrix(Vgt, i0, i1, j0, j1);

    free_f32matrix(Er_nir, i0, i1, j0, j1);
    free_f32matrix(Er_lin, i0, i1, j0, j1);
    free_f32matrix(Er_cub, i0, i1, j0, j1);

    free_ui8matrix(Ir_nir, i0, i1, j0, j1);
    free_ui8matrix(Ir_lin, i0, i1, j0, j1);
    free_ui8matrix(Ir_cub, i0, i1, j0, j1);

    return 0;
}
// ------------------------------------------
int HS_test_mono2_F32(int argc, char *argv[])
// ------------------------------------------
{
    int gt = 1; // 1 = synthetic, 0 = coadtguard
    int quantif = 0;

    int height, width;  // image size
    int i0, i1, j0, j1; // image intervals
    int b = 4;          // border

    float32 q; // quantification factor
    float32 alpha2;

    double ae, ee;                 // angle and endpoint errors
    double ie_nir, ie_lin, ie_cub; // interpolation error for {nearest, linear, cubic} interpolation

    char *src_path = "/nfs/home/sasl/eleves/main/3800549/Bureau/OpticalFlow_HS/OpticalFlow_HS/config/";
    char *dst_path = "/nfs/home/sasl/eleves/main/3800549/Bureau/OpticalFlow_HS/OpticalFlow_HS/config/";

    char *filename0;
    char *filename1;

    char complete_filename [1024];
    char complete_filename0[1024];
    char complete_filename1[1024];

    uint8 **I0;
    uint8 **I1;

    float32 **E0, **E1;       // pair of images (t), (t+1)
    float32 **Ex, **Ey, **Et; // {x, y, t} gradient

    float32 **U0, **V0; // input flow
    float32 **U1, **V1; // output flow
    float32 **Um, **Vm; // average flow (moyen)

    float32 **Ugt, **Vgt; // ground-true

    float32 **Er_nir, **Er_lin, **Er_cub;
    uint8   **Ir_nir, **Ir_lin, **Ir_cub;

    puts("-------------------------");
    puts("--- HS_test_mono2_F32 ---");
    puts("-------------------------");

    if(gt) {
        //synthetic cosine/sine
        filename0 = "I0.pgm";
        filename1 = "I1.pgm";
    } else {
        filename0 = "coastguard060.pgm";
        filename1 = "coastguard061.pgm";
        filename1 = "coastguard060q.pgm"; // quart de mvt de synthese
    }

    generate_path_filename(src_path, filename0, complete_filename0);
    generate_path_filename(src_path, filename1, complete_filename1);

    // load & alloc
    I0 = LoadPGM_ui8matrix(complete_filename0, &i0, &i1, &j0, &j1);
    I1 = LoadPGM_ui8matrix(complete_filename1, &i0, &i1, &j0, &j1);

    generate_path_filename(dst_path, "verif0.pgm", complete_filename0);
    generate_path_filename(dst_path, "verif1.pgm", complete_filename1);
    SavePGM_ui8matrix(I0, i0, i1, j0, j1, complete_filename0);
    SavePGM_ui8matrix(I1, i0, i1, j0, j1, complete_filename1);

    height = i1 - i0 + 1;
    width  = j1 - j0 + 1;
    printf("size = %d x %d\n", width, height);

    // alloc
    E0 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    E1 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Ex = f32matrix(i0, i1, j0, j1);
    Ey = f32matrix(i0, i1, j0, j1);
    Et = f32matrix(i0, i1, j0, j1);

    U0 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    V0 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    U1 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    V1 = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Um = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    Vm = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Ugt = f32matrix(i0, i1, j0, j1);
    Vgt = f32matrix(i0, i1, j0, j1);

    Er_nir = f32matrix(i0, i1, j0, j1);
    Er_lin = f32matrix(i0, i1, j0, j1);
    Er_cub = f32matrix(i0, i1, j0, j1);

    Ir_nir = ui8matrix(i0, i1, j0, j1);
    Ir_lin = ui8matrix(i0, i1, j0, j1);
    Ir_cub = ui8matrix(i0, i1, j0, j1);

    // init
    zero_f32matrix(Ex, i0, i1, j0, j1);
    zero_f32matrix(Ey, i0, i1, j0, j1);
    zero_f32matrix(Et, i0, i1, j0, j1);

    zero_f32matrix(E0, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_f32matrix(E1, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_f32matrix(U0, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_f32matrix(V0, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_f32matrix(U1, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_f32matrix(V1, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_f32matrix(Um, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_f32matrix(Vm, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_f32matrix(Ugt, i0, i1, j0, j1);
    zero_f32matrix(Vgt, i0, i1, j0, j1);

    if(gt) {
        set_f32matrix(Ugt, i0, i1, j0, j1, 0.75f);
        set_f32matrix(Vgt, i0, i1, j0, j1, 0.25f);
    }

    alpha2 = 0.635f; // mauvaise convergence
    alpha2 = 0.05f;  // bonne convergence
    //alpha2 = 0.04f;

    // convert - quantification 8 bits
    q = 255.0f; // attention: pas 256

    if(quantif) {
        // soit on normalise/quantifie les images sur [0,1]
        convert_ui8matrix_f32matrix_quantif(I0, i0, i1, j0, j1, q, E0);
        convert_ui8matrix_f32matrix_quantif(I1, i0, i1, j0, j1, q, E1);
    } else {
        // soit on re-quantifie alpha2
        convert_ui8matrix_f32matrix_quantif(I0, i0, i1, j0, j1, 1.0f, E0);
        convert_ui8matrix_f32matrix_quantif(I1, i0, i1, j0, j1, 1.0f, E1);
        alpha2 *= q*q;
    }

    // calc step #0
    hs_step0_F32(E0, E1, height, width, alpha2, Ex, Ey, Et, U1, V1);

    if(gt) {
        display_f32matrix(E0, height/4, 3*height/4, width/4, 3*width/4, "%8.2f", "E0");
        display_f32matrix(E1, height/4, 3*height/4, width/4, 3*width/4, "%8.2f", "E1");

        display_f32matrix(Ex, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "Ex");
        display_f32matrix(Ey, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "Ey");
        display_f32matrix(Et, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "Et");
    }

    ae = 0.0; ee = 0.0;
    for(int iter = 0; iter < 500; iter++) {

        hs_1step_F32(U0, V0, Ex, Ey, Et, height, width, alpha2, Um, Vm, U1, V1);
        if(gt) {
            ae = of_calcAngularError_F32 (U0, V0, Ugt, Vgt, height, width);
            ee = of_calcEndpointError_F32(U0, V0, Ugt, Vgt, height, width);
            printf("iter %3d angle = %8.3f  norm = %8.3f \n", iter, ae, ee);
        }
        SWAP_F32(U0, U1);
        SWAP_F32(V0, V1);
    }
    if(gt) {
        display_f32matrix(U0, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "U");
        display_f32matrix(V0, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "V");
    }

    generate_path_filename(dst_path, "U.txt", complete_filename0); write_f32matrix(U0, i0, i1, j0, j1, "%8.3f", complete_filename0);
    generate_path_filename(dst_path, "V.txt", complete_filename1); write_f32matrix(V0, i0, i1, j0, j1, "%8.3f", complete_filename1);

    interpNearest_F32(E1, U0, V0, height, width, Er_nir);
    interpLinear_F32 (E1, U0, V0, height, width, Er_lin);
    interpCubic_F32  (E1, U0, V0, height, width, Er_cub);

    if(quantif) {
        convert_f32matrix_ui8matrix_quantif(Er_nir, i0, i1, j0, j1, q, Ir_nir); // attention: pas 256
        convert_f32matrix_ui8matrix_quantif(Er_lin, i0, i1, j0, j1, q, Ir_lin);
        convert_f32matrix_ui8matrix_quantif(Er_cub, i0, i1, j0, j1, q, Ir_cub);
    } else {
        // impossible d'appeler directement la fonction convert_quantif, car clamping a q
        //convert_f32matrix_ui8matrix_quantif(Er_nir, i0, i1, j0, j1, 1.0f, Ir_nir);
        //convert_f32matrix_ui8matrix_quantif(Er_lin, i0, i1, j0, j1, 1.0f, Ir_lin);
        //convert_f32matrix_ui8matrix_quantif(Er_cub, i0, i1, j0, j1, 1.0f, Ir_cub);

        convert_f32matrix_ui8matrix_sat(Er_nir, i0, i1, j0, j1, Ir_nir);
        convert_f32matrix_ui8matrix_sat(Er_lin, i0, i1, j0, j1, Ir_lin);
        convert_f32matrix_ui8matrix_sat(Er_cub, i0, i1, j0, j1, Ir_cub);
    }

    display_f32matrix(Er_cub, height/4, 3*height/4, width/4, 3*width/4, "%6.3f ", "Er_cub");
    display_ui8matrix(Ir_cub, height/4, 3*height/4, width/4, 3*width/4, "%6d ",   "Ir_cub");

    generate_path_filename(dst_path, "Irec_nir.pgm", complete_filename); SavePGM_ui8matrix(Ir_nir, i0, i1, j0, j1, complete_filename);
    generate_path_filename(dst_path, "Irec_lin.pgm", complete_filename); SavePGM_ui8matrix(Ir_lin, i0, i1, j0, j1, complete_filename);
    generate_path_filename(dst_path, "Irec_cub.pgm", complete_filename); SavePGM_ui8matrix(Ir_cub, i0, i1, j0, j1, complete_filename);

    generate_path_filename(dst_path, "Irec_nir.txt", complete_filename); write_ui8matrix_number(Ir_nir, i0, i1, j0, j1, "%4d", complete_filename);
    generate_path_filename(dst_path, "Irec_lin.txt", complete_filename); write_ui8matrix_number(Ir_lin, i0, i1, j0, j1, "%4d", complete_filename);
    generate_path_filename(dst_path, "Irec_cub.txt", complete_filename); write_ui8matrix_number(Ir_cub, i0, i1, j0, j1, "%4d", complete_filename);

    if(gt) {
        ae = of_calcAngularError_F32 (U0, V0, Ugt, Vgt, height, width);
        ee = of_calcEndpointError_F32(U0, V0, Ugt, Vgt, height, width);
        printf("error angle(deg) = %.3f  error L2norm(pixel) = %.3f\n", ae, ee);
    }
    ie_nir = of_calcInterpolationError_F32(E0, Er_nir, height, width);
    ie_lin = of_calcInterpolationError_F32(E0, Er_lin, height, width);
    ie_cub = of_calcInterpolationError_F32(E0, Er_cub, height, width);

    if(quantif) {
        printf("interp_nir(gray) = %.3f  interp_lin(gray) = %.3f  interp_cub(gray) = %.3f\n", q*ie_nir, q*ie_lin, q*ie_cub);
    } else {
        printf("interp_nir(gray) = %.3f  interp_lin(gray) = %.3f  interp_cub(gray) = %.3f\n", ie_nir, ie_lin, ie_cub);
    }

    // free
    puts("free");
    free_ui8matrix(I0, i0, i1, j0, j1);
    free_ui8matrix(I1, i0, i1, j0, j1);

    free_f32matrix(E0, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(E1, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Ex, i0, i1, j0, j1);
    free_f32matrix(Ey, i0, i1, j0, j1);
    free_f32matrix(Et, i0, i1, j0, j1);

    free_f32matrix(U0, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(V0, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(U1, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(V1, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Um, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(Vm, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Ugt, i0, i1, j0, j1);
    free_f32matrix(Vgt, i0, i1, j0, j1);

    free_f32matrix(Er_nir, i0, i1, j0, j1);
    free_f32matrix(Er_lin, i0, i1, j0, j1);
    free_f32matrix(Er_cub, i0, i1, j0, j1);

    free_ui8matrix(Ir_nir, i0, i1, j0, j1);
    free_ui8matrix(Ir_lin, i0, i1, j0, j1);
    free_ui8matrix(Ir_cub, i0, i1, j0, j1);

    return 0;
}
// -----------------------------------------
int HS_test_mono1_FQ(int argc, char *argv[])
// -----------------------------------------
{
#ifdef ENABLE_FIXED_POINT
    int height, width;  // image size
    int i0, i1, j0, j1; // image intervals
    int b = 4;          // border

    int q_image, q_derivative, q_speed, q_interp; // speed quantification
    sint64 Q_image, Q_derivative, Q_speed, Q_interp;

    float32 alpha2; // valeur F32 pour des images de dynamique 1: [0..1] (et pas sur 0..255)
    float32 alpha2_f; // pour des images F32 quantifiees
    sint32  alpha2_q; // pour des images I16 quantifiees

    //double ie;     // interpolation error for cubic interpolation

    char *src_path = "/nfs/home/sasl/eleves/main/3800549/Bureau/OpticalFlow_HS/OpticalFlow_HS/config/";
    char *dst_path = "/Users/lacas/Code/OpticalFlow/result/";

    char *filename0;
    char *filename1;

    char complete_filename [1024];
    char complete_filename0[1024];
    char complete_filename1[1024];

    // 8-bit input
    uint8 **I0;
    uint8 **I1;

    uint16 **E0, **E1;       // pair of images (t), (t+1)
    sint16 **Ex, **Ey, **Et; // {x, y, t} gradient

    sint32 **U0, **V0; // input flow
    sint32 **U1, **V1; // output flow
    sint32 **Um, **Vm; // average flow (moyen)

    uint16 **Er;

    float32 **E0_f, **E1_f;        // pair of images (t), (t+1)
    float32 **Ex_f, **Ey_f, **Et_f; // {x, y, t} gradient

    float32 **U0_f, **V0_f; // input flow
    float32 **U1_f, **V1_f; // output flow
    float32 **Um_f, **Vm_f; // average flow (moyen)

    float32 **cU, **cV;

    float32 **Er_f;

    // 8-bit output
    uint8 **Ir, **Ir_f;

    //double eeUV; // end-point error entre 2 vecteurs
    //double mse_U, mse_V;
    double mse_FF, mse_IF, mse_II; // erreur d'interpolation cubique
    double mse_FL, mse_IL; // erreur d'interpolation lineaire
    //double er_U, er_V;

    puts("-------------------------");
    puts("--- HS_test_mono1_F+Q ---");
    puts("-------------------------");

    // chargement 8 bits (prevoir 16 pour la suite)

    filename0 = "I0.pgm"; filename1 = "I1.pgm"; //synthetic cosine/sine
    filename0 = "coastguard060.pgm"; filename1 = "coastguard060q.pgm"; // quart de mouvement de synthese

    generate_path_filename(src_path, filename0, complete_filename0);
    generate_path_filename(src_path, filename1, complete_filename1);

    // load & alloc
    I0 = LoadPGM_ui8matrix(complete_filename0, &i0, &i1, &j0, &j1);
    I1 = LoadPGM_ui8matrix(complete_filename1, &i0, &i1, &j0, &j1);

    height = i1 - i0 + 1;
    width  = j1 - j0 + 1;

    int i0d, i1d, j0d, j1d;

    i0d = 1*height/4; i1d = 3*height/4; j0d = 1*width/4; j1d = 1*width/4; // 1/4 .. 3/4
    i0d = 1*height/8; i1d = 2*height/8; j0d = 1*width/8; j1d = 2*width/8; // 1/8 .. 2/8
    i0d = 10; i1d = i0d + 7; j0d = 10; j1d = j0d + 7; //

    Ir   = ui8matrix(i0, i1, j0, j1);
    Ir_f = ui8matrix(i0, i1, j0, j1);

    printf("size = %d x %d\n", width, height);

    // alloc INT
    E0 = ui16matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    E1 = ui16matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Ex = si16matrix(i0, i1, j0, j1);
    Ey = si16matrix(i0, i1, j0, j1);
    Et = si16matrix(i0, i1, j0, j1);

    U0 = si32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    V0 = si32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    U1 = si32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    V1 = si32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Um = si32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    Vm = si32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Er = ui16matrix(i0, i1, j0, j1);

    // alloc FLOAT
    E0_f = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    E1_f = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Ex_f = f32matrix(i0, i1, j0, j1);
    Ey_f = f32matrix(i0, i1, j0, j1);
    Et_f = f32matrix(i0, i1, j0, j1);

    U0_f = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    V0_f = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    U1_f = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    V1_f = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Um_f = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    Vm_f = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    Er_f = f32matrix(i0, i1, j0, j1);

    // float for comparison
    cU = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);
    cV = f32matrix(i0 - b, i1 + b, j0 - b, j1 + b);

    // zero

    zero_ui16matrix(E0, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_ui16matrix(E1, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_si16matrix(Ex, i0, i1, j0, j1);
    zero_si16matrix(Ey, i0, i1, j0, j1);
    zero_si16matrix(Et, i0, i1, j0, j1);

    zero_si32matrix(U0, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_si32matrix(V0, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_si32matrix(U1, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_si32matrix(V1, i0 - b, i1 + b, j0 - b, j1 + b);

    zero_si32matrix(Um, i0 - b, i1 + b, j0 - b, j1 + b);
    zero_si32matrix(Vm, i0 - b, i1 + b, j0 - b, j1 + b);

    // conversion pour quantification
    convert_ui8matrix_ui16matrix(I0, i0, i1, j0, j1, E0);
    convert_ui8matrix_ui16matrix(I1, i0, i1, j0, j1, E1);

    convert_ui8matrix_f32matrix(I0, i0, i1, j0, j1, E0_f);  // attention: pas 256
    convert_ui8matrix_f32matrix(I1, i0, i1, j0, j1, E1_f);

    //display_f32matrix(E0, height/4, 3*height/4, width/4, 3*width/4, "%8.2f", "E0");
    //display_f32matrix(E1, height/4, 3*height/4, width/4, 3*width/4, "%8.2f", "E1");

    // calc

    // quantification 'basique'
    q_image      = 8;
    q_derivative = 8;
    q_speed      = 8;
    q_interp     = 8;

    // quantification 'embedded'
    /*q_image      = 8;
    q_derivative = 7;
    q_speed      = 4;
    q_interp     = 5;*/

    // base  2
    //q = 2; // base 10

    Q_image      = CALC_Q(q_image);
    Q_derivative = CALC_Q(q_derivative);
    Q_speed      = CALC_Q(q_speed);
    Q_interp     = CALC_Q(q_interp);

    alpha2 = 0.635f; // mauvaise convergence
    alpha2 = 0.05f;  // bonne convergence
    alpha2 = 0.04f;
    alpha2 = 0.02f;
    alpha2 = 0.01f;
    //alpha2 = 0.005f;
    //alpha2 = 0.002f;
    //alpha2 = 0.001f;

    alpha2_f = alpha2 * Q_image * Q_image;
    alpha2_q = (int) roundf(alpha2 * Q_derivative * Q_derivative);

    printf("qi = %d Qi = %lld\n", q_image,      Q_image);
    printf("qd = %d Qd = %lld\n", q_derivative, Q_derivative);
    printf("qs = %d Qs = %lld\n", q_speed,      Q_speed);

    // ================ //
    // === prologue === //
    // ================ //

    hs_step0_F32 (E0_f, E1_f, height, width, alpha2_f, Ex_f, Ey_f, Et_f, U0_f, V0_f);
    //hs_step0_Qstd(E0, E1, q_image, height, width, alpha2_q, Ex, Ey, Et, q_derivative, U0, V0, q_speed);
    hs_step0_FQstd(E0, E1, q_image, height, width, alpha2_q, Ex, Ey, Et, q_derivative, U0, V0, q_speed);

    //hs_step0_Qmax(E0, E1, q_image, height, width, alpha2_q, Ex, Ey, Et, q_derivative, U0, V0, q_speed);

    puts("AVANT");

    display_ui16matrix(E0, i0d, i1d, j0d, j1d, "%8d", "E0 I16");
    display_ui16matrix(E1, i0d, i1d, j0d, j1d, "%8d", "E1 I16");

    display_f32matrix(E0_f, i0d, i1d, j0d, j1d, "%8.0f", "E0 F32");
    display_f32matrix(E1_f, i0d, i1d, j0d, j1d, "%8.0f", "E1 F32");

    display_si16matrix(Ex, i0d, i1d, j0d, j1d, "%8d", "Ex I16");
    display_si16matrix(Ey, i0d, i1d, j0d, j1d, "%8d", "Ey I16");
    display_si16matrix(Et, i0d, i1d, j0d, j1d, "%8d", "Et I16");

    // Qmax(derivative) = 4 * Qstd(derivative)
    display_si16matrix_quantif(Ex, i0d, i1d, j0d, j1d, 4, "%8.3f", "Ex Q16");
    display_si16matrix_quantif(Ey, i0d, i1d, j0d, j1d, 4, "%8.3f", "Ey Q16");
    display_si16matrix_quantif(Et, i0d, i1d, j0d, j1d, 4, "%8.3f", "Et Q16");/**/

    display_f32matrix(Ex_f, i0d, i1d, j0d, j1d, "%8.3f", "Ex F32");
    display_f32matrix(Ey_f, i0d, i1d, j0d, j1d, "%8.3f", "Ey F32");
    display_f32matrix(Et_f, i0d, i1d, j0d, j1d, "%8.3f", "Et F32");

    display_si32matrix_quantif(U0, i0d, i1d, j0d, j1d, Q_speed, "%8.3f", "U0 I32");
    display_si32matrix_quantif(V0, i0d, i1d, j0d, j1d, Q_speed, "%8.3f", "V0 I32");

    display_f32matrix(U0_f, i0d, i1d, j0d, j1d, "%8.3f", "U0 F32");
    display_f32matrix(V0_f, i0d, i1d, j0d, j1d, "%8.3f", "V0 F32");/**/

    generate_path_filename(dst_path, "U_.txt", complete_filename);  write_si32matrix(U0, i0, i1, j0, j1, "%6d", complete_filename);
    generate_path_filename(dst_path, "V_.txt", complete_filename);  write_si32matrix(V0, i0, i1, j0, j1, "%6d", complete_filename);

    //return 0;

    // ------------------
    // -- erreur initiale
    // ------------------
    mse_FF = of_calcMeanSquareError_F32(E0_f, E1_f, height, width);
    printf("iter = %02d mse_FF = %8.3f\n\n", 0, mse_FF);

    // ----------------------------------
    // -- erreur apres iteration initiale
    // ----------------------------------
    printf("iter = %3d ", 0);

    convert_si32matrix_f32matrix_quantif(U0, 0, height-1, 0, width-1, Q_speed, cU);
    convert_si32matrix_f32matrix_quantif(V0, 0, height-1, 0, width-1, Q_speed, cV);

    // difference des vecteurs vitesses
    /*eeUV = of_calcEndpointError_F32(cU, cV, U0_f, V0_f, height, width);
    printf("ee = %8.3f   ", eeUV);*/

    /*mse_U = of_calcMeanSquareError_F32(cU, U0_f, height, width);
    mse_V = of_calcMeanSquareError_F32(cV, V0_f, height, width);
    printf("mse_U = %8.3f   mse_V = %8.3f   ", mse_U, mse_V);*/

    // calcul F + interpolation cubique F
    interpCubic_F32(E1_f, U0_f, V0_f, height, width, Er_f);
    mse_FF  = of_calcMeanSquareError_F32(E0_f, Er_f, height, width);

    // calcul I + interpolation cubique F(I)
    interpCubic_F32(E1_f, cU, cV, height, width, Er_f);
    mse_IF = of_calcMeanSquareError_F32(E0_f, Er_f, height, width);

    // calcul I + interpolation cubique I
    interpCubic_I16_I32(E1, q_image, U0, V0, q_speed, height, width, Er);
    convert_ui16matrix_f32matrix(Er, 0, height-1, 0, width-1, Er_f);
    mse_II = of_calcMeanSquareError_F32(Er_f, E0_f, height, width);


    // calcul F + interpolation lineaire F
    interpLinear_F32(E1_f, U0_f, V0_f, height, width, Er_f);
    mse_FL = of_calcMeanSquareError_F32(Er_f, E0_f, height, width);

    // calcul I + interpolation lineaire I
    interpLinear_I16_I32(E1, q_image, U0, V0, q_speed, height, width, Er);
    convert_ui16matrix_f32matrix(Er, 0, height-1, 0, width-1, Er_f);
    mse_IL = of_calcMeanSquareError_F32(Er_f, E0_f, height, width);

    //printf("iter = %2d mse_F = %8.3f   mse_I = %8.3f\n", 0, mse_F, mse_I);
    printf("mse_FF = %8.3f   mse_IF = %8.3f   mse_II = %8.3f\n", mse_FF, mse_IF, mse_II);
    //printf("mse_FF = %6.3f   mse_IF = %6.3f   mse_II = %6.3f  |  mse_FLI = %6.3f   mse_FI = %6.3f", mse_FF, mse_IF, mse_II, mse_FL, mse_IL);

    // erreur relative
    //er_U = og_calcRelativeError_F32(U0_f, cU, height, width);
    //er_V = og_calcRelativeError_F32(V0_f, cV, height, width);
    //printf("iter = %2d er_U = %8.3f   er_V = %8.3f\n", 0, er_U, er_V);

    putchar('\n');

    // ============== //
    // === boucle === //
    // ============== //

    for(int iter = 1; iter < 20; iter++) {

        printf("iter = %3d ", iter);

        hs_1step_F32(U0_f, V0_f, Ex_f, Ey_f, Et_f, height, width, alpha2_f, Um_f, Vm_f, U1_f, V1_f);
        //hs_1step_Qstd(U0, V0, q_speed, Ex, Ey, Et, q_derivative, height, width, alpha2_q, Um, Vm, U1, V1);
        hs_1step_FQstd(U0, V0, q_speed, Ex, Ey, Et, q_derivative, height, width, alpha2_q, Um, Vm, U1, V1);
        //hs_1step_Qmax(U0, V0, q_speed, Ex, Ey, Et, q_derivative, height, width, alpha2_q, Um, Vm, U1, V1);

        convert_si32matrix_f32matrix_quantif(U1, 0, height-1, 0, width-1, Q_speed, cU);
        convert_si32matrix_f32matrix_quantif(V1, 0, height-1, 0, width-1, Q_speed, cV);

        //er_U = og_calcRelativeError_F32(U1_f, cU, height, width);
        //er_V = og_calcRelativeError_F32(V1_f, cV, height, width);
        //printf("iter = %2d er_U = %8.3f   er_V = %8.3f\n", iter, er_U, er_V);

        // difference des vecteurs vitesses
        /*eeUV = of_calcEndpointError_F32(cU, cV, U1_f, V1_f, height, width);
        printf("ee = %8.3f   ", eeUV);*/

        /*mse_U = of_calcMeanSquareError_F32(cU, U1_f, height, width);
        mse_V = of_calcMeanSquareError_F32(cV, V1_f, height, width);
        printf("mse_U = %10.5f   mse_V = %10.5f   ", mse_U, mse_V);*/

        // calcul F + interpolation cubique F
        interpCubic_F32(E1_f, U1_f, V1_f, height, width, Er_f);
        mse_FF = of_calcMeanSquareError_F32(Er_f, E0_f, height, width);

        // calcul I + interpolation cubique F(I)
        interpCubic_F32(E1_f, cU, cV, height, width, Er_f);
        mse_IF = of_calcMeanSquareError_F32(Er_f, E0_f, height, width);

        // calcul I + interpolation cubique I
        interpCubic_I16_I32(E1, q_image, U1, V1, q_speed, height, width, Er);
        convert_ui16matrix_f32matrix(Er, 0, height-1, 0, width-1, Er_f);
        mse_II = of_calcMeanSquareError_F32(Er_f, E0_f, height, width);

        // calcul F + interpolation lineaire F
        interpLinear_F32(E1_f, U1_f, V1_f, height, width, Er_f);
        mse_FL = of_calcMeanSquareError_F32(Er_f, E0_f, height, width);

        // calcul I + interpolation lineaire I
        interpLinear_I16_I32(E1, q_image, U1, V1, q_speed, height, width, Er);
        convert_ui16matrix_f32matrix(Er, 0, height-1, 0, width-1, Er_f);
        mse_IL = of_calcMeanSquareError_F32(Er_f, E0_f, height, width);

        printf("mse_FF = %6.3f   mse_IF = %6.3f   mse_II = %6.3f", mse_FF, mse_IF, mse_II);
        //printf("mse_FF = %6.3f   mse_IF = %6.3f   mse_II = %6.3f  |  mse_FL = %6.3f   mse_IL = %6.3f", mse_FF, mse_IF, mse_II, mse_FL, mse_IL);

        SWAP_F32(U0_f, U1_f);
        SWAP_F32(V0_f, V1_f);

        SWAP_SI32(U0, U1);
        SWAP_SI32(V0, V1);

        putchar('\n');
    }

    // ================ //
    // === epilogue === //
    // ================ //

    puts("APRES");
    //display_si64matrix_quantif(U0, height/4, 3*height/4, width/4, 3*width/4, q, "%8.3f", "U0 I64");
    //display_si64matrix_quantif(V0, height/4, 3*height/4, width/4, 3*width/4, q, "%8.3f", "V0 I64");

    //display_f32matrix(U0_f, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "U0 F32");
    //display_f32matrix(V0_f, height/4, 3*height/4, width/4, 3*width/4, "%8.3f", "V0 F32");

    //display_si64matrix(U0, height/4, 3*height/4, width/4, 3*width/4, "%8d", "U");
    //display_si64matrix(V0, height/4, 3*height/4, width/4, 3*width/4, "%8d", "V");

    generate_path_filename(dst_path, "Ui.txt", complete_filename); write_si32matrix(U0, i0, i1, j0, j1, "%6d", complete_filename);
    generate_path_filename(dst_path, "Vi.txt", complete_filename); write_si32matrix(V0, i0, i1, j0, j1, "%6d", complete_filename);

    generate_path_filename(dst_path, "Uf.txt", complete_filename); write_f32matrix(U0_f, i0, i1, j0, j1, "%6.2f", complete_filename);
    generate_path_filename(dst_path, "Vf.txt", complete_filename); write_f32matrix(V0_f, i0, i1, j0, j1, "%6.2f", complete_filename);

    // interpolation virgule fixe + conversion sans quantif car input sur 8 bit pour le moment ...
    interpCubic_I16_I32(E1, q_image, U0, V0, q_speed, height, width, Er);
    convert_ui16matrix_ui8matrix_sat(Er, i0, i1, j0, j1, Ir);
    generate_path_filename(dst_path, "Irec.pgm", complete_filename); SavePGM_ui8matrix(Ir, i0, i1, j0, j1, complete_filename);

    // interpolation virgule flottante + conversion avec saturation
    interpCubic_F32(E1_f, U0_f, V0_f, height, width, Er_f);
    convert_f32matrix_ui8matrix_sat(Er_f, i0, i1, j0, j1, Ir_f);
    generate_path_filename(dst_path, "Frec.pgm", complete_filename); SavePGM_ui8matrix(Ir_f, i0, i1, j0, j1, complete_filename);

    // interpolation virgule flottante depuis (U,V) en virgule fixe (pour tester qualite interp virgule fixe)
    convert_si32matrix_f32matrix_quantif(U0, 0, height-1, 0, width-1, Q_speed, cU);
    convert_si32matrix_f32matrix_quantif(V0, 0, height-1, 0, width-1, Q_speed, cV);
    interpCubic_F32(E1_f, cU, cV, height, width, Er_f);
    convert_f32matrix_ui8matrix_sat(Er_f, i0, i1, j0, j1, Ir_f);
    generate_path_filename(dst_path, "IFrec.pgm", complete_filename); SavePGM_ui8matrix(Ir_f, i0, i1, j0, j1, complete_filename);

    //display_ui8matrix(Er, height/4, 3*height/4, width/4, 3*width/4, "%6d ", "Er");
    //display_ui8matrix(Ir, height/4, 3*height/4, width/4, 3*width/4, "%6d ", "Ir");

    //generate_path_filename(dst_path, "Erec.txt", complete_filename); write_ui8matrix_number(Er, i0, i1, j0, j1, "%4d", complete_filename);

    //ie = of_calcInterpolationError_I64(E0, Er, height, width);

    //printf("interp = %.3f\n", ie);

    printf("=> autant il est possible d'obtenir une interpolation cubique precise en Q (trois premieres colonnes)\n");
    printf("=> autant il est impossible d'obtenir une interpolation lineaire precise en Q (deux dernieres colonnes)\n");
    // free
    puts("free");
    free_ui8matrix(I0, i0, i1, j0, j1);
    free_ui8matrix(I1, i0, i1, j0, j1);

    free_ui8matrix(Ir, i0, i1, j0, j1);
    free_ui8matrix(Ir_f, i0, i1, j0, j1);

    // Int
    free_ui16matrix(E0, i0 - b, i1 + b, j0 - b, j1 + b);
    free_ui16matrix(E1, i0 - b, i1 + b, j0 - b, j1 + b);

    free_si16matrix(Ex, i0, i1, j0, j1);
    free_si16matrix(Ey, i0, i1, j0, j1);
    free_si16matrix(Et, i0, i1, j0, j1);

    free_si32matrix(U0, i0 - b, i1 + b, j0 - b, j1 + b);
    free_si32matrix(V0, i0 - b, i1 + b, j0 - b, j1 + b);

    free_si32matrix(U1, i0 - b, i1 + b, j0 - b, j1 + b);
    free_si32matrix(V1, i0 - b, i1 + b, j0 - b, j1 + b);

    free_si32matrix(Um, i0 - b, i1 + b, j0 - b, j1 + b);
    free_si32matrix(Vm, i0 - b, i1 + b, j0 - b, j1 + b);

    free_ui16matrix(Er, i0, i1, j0, j1);

    // Float
    free_f32matrix(E0_f, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(E1_f, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Ex_f, i0, i1, j0, j1);
    free_f32matrix(Ey_f, i0, i1, j0, j1);
    free_f32matrix(Et_f, i0, i1, j0, j1);

    free_f32matrix(U0_f, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(V0_f, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(U1_f, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(V1_f, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Um_f, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(Vm_f, i0 - b, i1 + b, j0 - b, j1 + b);

    free_f32matrix(Er_f, i0, i1, j0, j1);

    free_f32matrix(cU, i0 - b, i1 + b, j0 - b, j1 + b);
    free_f32matrix(cV, i0 - b, i1 + b, j0 - b, j1 + b);
#endif // ENABLE_FIXED_POINT
    return 0;
}
// --------------------------------------------------
int HS_test_hierarchique1_F32(int argc, char *argv[])
// --------------------------------------------------
{
    int h, w, height, width;  // image size
    int i0, i1, j0, j1; // image intervals
    int nb_digit = 2;
    int iter,  nb_iter  = 200;
    int level, nb_level = 3;
    int border = 3; // for interpolation

    float32 alpha2 = 0.005f;

    double mse;

    char *src_path = "/nfs/home/sasl/eleves/main/3800549/Bureau/OpticalFlow_HS/OpticalFlow_HS/config/";
    char *dst_path = "/Users/lacas/Code/OpticalFlow/resultHierarchique/";

    char *filename0;
    char *filename1;

    char complete_filename [1024];
    char complete_filename0[1024];
    char complete_filename1[1024];

    uint8   **I0, **I1, **I1r;
    float32 **K; // Burt Kernel

    float32 ***F0, ***F1, ***F1r, ***F1t;         // pair of images (t), (t+1), reconstruction, reconstruction temporaire
    float32 ***Fx, ***Fy, ***Ft;                  // {x, y, t} gradient
    float32 ***dU, ***dV;
    float32 ***Um, ***Vm;
    float32 ***U, ***V;
    float32 ***Ut, ***Vt;
    float32 **W;

    puts("---------------------------------");
    puts("--- HS_test_hierarchique1_F32 ---");
    puts("---------------------------------");

    // pour la paire d'images coastguard060 + coastguard061, 3 niveaux, 100 iters
    //alpha2 = 0.100f; // mse(100)=52.81
    //alpha2 = 0.050f; // mse(100)=43.54
    //alpha2 = 0.020f; // mse(100)=36.42
    //alpha2 = 0.010f; // mse(100)=33.95
    //alpha2 = 0.008f; // mse(100)=33.58
    //alpha2 = 0.007f; // mse(100)=33.46
    //alpha2 = 0.006f; // mse(100)=33.43
    //alpha2 = 0.005f; // mse(100)=33.53
    //alpha2 = 0.004f; // mse(100)=33.86
    //alpha2 = 0.003f; // mse(100)=35.42
    //alpha2 = 0.002f; // mse(100)=38.24
    //alpha2 = 0.001f; // mse(100)=43.77

    // pour la paire d'images coastguard150 + coastguard151, 3 niveaux, 100 iters

    //alpha2 = 0.100f; // mse = 46.64
    //alpha2 = 0.050f; // mse = 35.87
    //alpha2 = 0.020f; // mse = 25.14
    //alpha2 = 0.010f; // mse = 20.31
    //alpha2 = 0.008f; // mse = 19.34
    //alpha2 = 0.007f; // mse = 18.87
    //alpha2 = 0.006f; // mse = 18.43 <- best de (060,061)
    //alpha2 = 0.005f; // mse = 18.06
    //alpha2 = 0.004f; // mse = 17.81
    //alpha2 = 0.003f; // mse = 17.80 <- best
    //alpha2 = 0.002f; // mse = 18.31
    //alpha2 = 0.001f; // mse = 20.47

    alpha2 = 0.0100f;
    alpha2 = 0.0080f;
    //alpha2 = 0.0060f;
    alpha2 = 0.0050f;
    alpha2 = 0.0040f;
    //alpha2 = 0.0030f;
    //alpha2 = 0.0020f;
    //alpha2 = 0.0010f;
    //alpha2 = 0.0005f;
    alpha2 = 0.0009f;
    alpha2 = 0.0008f;
    alpha2 = 0.0007f;
    alpha2 = 0.0001f;
    alpha2 = 0.0500f;

    nb_level = 3;
    nb_iter  = 20;

    //filename0 = "I0.pgm"; filename1 = "I1.pgm"; //synthetic cosine/sine
    //filename0 = "c60.pgm"; filename1 = "c60t.pgm";
    //filename0 = "meteor1_0495.pgm"; filename1 = "meteor1_0496.pgm";

    //filename0 = "coastguard060.pgm"; filename1 = "coastguard060t.pgm";
    //filename0 = "coastguard060.pgm"; filename1 = "coastguard061.pgm";

    //filename0 = "coastguard060.pgm"; filename1 = "coastguard060t.pgm";
    //filename0 = "coastguard060.pgm"; filename1 = "coastguard060q.pgm";

    //filename0 = "coastguard149.pgm"; filename1 = "coastguard150.pgm";
    filename0 = "coastguard150.pgm"; filename1 = "coastguard151.pgm"; // descend plus bas

    generate_path_filename(src_path, filename0, complete_filename0);
    generate_path_filename(src_path, filename1, complete_filename1);

    // load & alloc
    I0 = LoadPGM_ui8matrix(complete_filename0, &i0, &i1, &j0, &j1);
    I1 = LoadPGM_ui8matrix(complete_filename1, &i0, &i1, &j0, &j1);

    // full size
    height = i1 - i0 + 1;
    width  = j1 - j0 + 1;

    printf("size = %d x %d\n", width, height);
    printf("%s vs %s\n", filename0, filename1);
    printf("nb_level = %d\n", nb_level);
    printf("nb_iter = %d\n", nb_iter);
    printf("alpha = %.4f\n", alpha2);

    // -----------
    // -- alloc --
    // -----------

    I1r  = ui8matrix(0, height-1, 0, width-1);

    K = f32matrix(-2, +2, -2, +2); // 5x5 kernel for Burt Pyramid

    F0  = f32pyramid(nb_level, height, width, border);
    F1  = f32pyramid(nb_level, height, width, border);
    F1r = f32pyramid(nb_level, height, width, border);
    F1t = f32pyramid(nb_level, height, width, border);

    Fx = f32pyramid(nb_level, height, width, border);
    Fy = f32pyramid(nb_level, height, width, border);
    Ft = f32pyramid(nb_level, height, width, border);

    dU = f32pyramid(nb_level, height, width, border);
    dV = f32pyramid(nb_level, height, width, border);

    Um = f32pyramid(nb_level, height, width, border);
    Vm = f32pyramid(nb_level, height, width, border);

    U = f32pyramid(nb_level, height, width, border);
    V = f32pyramid(nb_level, height, width, border);

    Ut = f32pyramid(nb_level, height, width, 0);
    Vt = f32pyramid(nb_level, height, width, 0);

    W = f32matrix(0, height-1, 0, width-1);

    // ----------
    // -- init --
    // ----------

    zero_f32pyramid(F0,  nb_level, height, width, border);
    zero_f32pyramid(F1,  nb_level, height, width, border);
    zero_f32pyramid(F1r, nb_level, height, width, border);
    zero_f32pyramid(F1t, nb_level, height, width, border);

    zero_f32pyramid(Fx, nb_level, height, width, border);
    zero_f32pyramid(Fy, nb_level, height, width, border);
    zero_f32pyramid(Ft, nb_level, height, width, border);

    zero_f32pyramid(dU, nb_level, height, width, border);
    zero_f32pyramid(dV, nb_level, height, width, border);

    zero_f32pyramid(Um, nb_level, height, width, border);
    zero_f32pyramid(Vm, nb_level, height, width, border);

    zero_f32pyramid(U, nb_level, height, width, border);
    zero_f32pyramid(V, nb_level, height, width, border);

    zero_f32pyramid(Ut, nb_level, height, width, 0);
    zero_f32pyramid(Vt, nb_level, height, width, 0);

    zero_f32matrix(W, 0, height-1, 0, width-1);

    // convert avec division par 255: [0:255] -> [0:1]
    //convert_ui8matrix_f32matrix_quantif(I0, i0, i1, j0, j1, 255, E0);  // attention: pas 256
    //convert_ui8matrix_f32matrix_quantif(I1, i0, i1, j0, j1, 255, E1);

    // -------------------------------- //
    // -- construction des pyramides -- //
    // -------------------------------- //

    alpha2 *= (255.f * 255.f);

    //for(level = 0; level < nb_level; level++) { Alpha2[level] = alpha2; }

    initBurtKernel_F32(0.4f, K);

    // prologue d'init

    level = 0;
    h = height >> level;
    w = width  >> level;

    convert_ui8matrix_f32matrix(I0, 0, h-1, 0, w-1, F0[level]);
    convert_ui8matrix_f32matrix(I1, 0, h-1, 0, w-1, F1[level]);

    // pour debug
    //convert_f32matrix_ui8matrix(F0[level], 0, h-1, 0, w-1, I0);
    //convert_f32matrix_ui8matrix(F1[level], 0, h-1, 0, w-1, I1);
    //generate_path_filename_k_ndigit_extension(dst_path, "I0_", level, nb_digit, "pgm", complete_filename0);
    //generate_path_filename_k_ndigit_extension(dst_path, "I1_", level, nb_digit, "pgm", complete_filename1);
    //SavePGM_ui8matrix(I0, 0, h-1, 0, w-1, complete_filename0);
    //SavePGM_ui8matrix(I1, 0, h-1, 0, w-1, complete_filename1);

    // -------------------
    // -- boucle d'init --
    // -------------------

    for(level = 1; level <= nb_level - 1 ; level++) {

        h = height >> level;
        w = width  >> level;

        // Burt pyramid
        convolution5x5_decimation_F32_F32_F32(F0[level-1], 2*h, 2*w, K, F0[level]);
        convolution5x5_decimation_F32_F32_F32(F1[level-1], 2*h, 2*w, K, F1[level]);

        // pour debug
        /*convert_f32matrix_ui8matrix(F0[level], 0, h-1, 0, w-1, I0);
        convert_f32matrix_ui8matrix(F1[level], 0, h-1, 0, w-1, I1);
        generate_path_filename_k_ndigit_extension(dst_path, "I0_", level, nb_digit, "pgm", complete_filename0);
        generate_path_filename_k_ndigit_extension(dst_path, "I1_", level, nb_digit, "pgm", complete_filename1);
        SavePGM_ui8matrix(I0, 0, h-1, 0, w-1, complete_filename0);
        SavePGM_ui8matrix(I1, 0, h-1, 0, w-1, complete_filename1);/**/

    } // level

    // ---------------- //
    // -- traitement -- //
    // ---------------- //

    // erreur initiale (full res)
    mse = of_calcMeanSquareError_F32(F0[0], F1[0], h, w);
    printf("level = %d iter = %3d mse = %8.3f\n", -1, -1, mse);

    // =====================================
    // == boucle pyramidale sur le niveau ==
    // =====================================

    for(level = nb_level - 1; level >= 0; level--) {

        iter = -1;
        printf("-----------------------------------\n");
        printf("level = %d\n", level);
        h = height >> level;
        w = width  >> level;
        i0 = 0; i1 = h-1;
        j0 = 0; j1 = w-1;

        // erreur initiale au niveau courant, avant le recallage
        mse = of_calcMeanSquareError_F32(F0[level], F1[level], h, w);
        printf("level = %d iter = %3d mse = %8.3f avant warp initial\n", level, iter, mse);

        // inutile normalement car write before read
        //zero_f32matrix(dU[level], 0-border, h-1+border, 0-border, w-1+border);
        //zero_f32matrix(dV[level], 0-border, h-1+border, 0-border, w-1+border);

        // ----------------------------------
        // -- recalage de debut de niveau --
        // ----------------------------------

        // recalage / compensation du mouvement avec les vecteurs vitesses calcules et upsamples du niveau precedent

        extendDataD2_F32(F1[level], h, w);
        //if(level < nb_level-1) { interpCubic_F32(F1[level], U[level], V[level], h, w, F1r[level]); }

        interpCubic_F32(F1[level], U[level], V[level], h, w, F1r[level]);

        mse = of_calcMeanSquareError_F32(F0[level], F1r[level], h, w);
        printf("level = %d iter = %3d mse = %8.3f apres warp initial\n", level, iter, mse);

        // iteration initiale avec (U,V)=0 => (Um,Vm)=0
        /*if(level == nb_level-1) {
            hs_step0_F32(F0[level], F1[level], h, w, alpha2, Fx[level], Fy[level], Ft[level], dU[level], dV[level]);
        }*/

        // ----------------------------------------------------------------------------------
        // -- step0: calcul des gradients (Ix,Iy,It) et init des vecteurs vitesses (dU,dV) --
        // ----------------------------------------------------------------------------------

        // naivement (faux)
        //put("calcul avec reset de (U,V) a chaque niveau et initialisation via derivees");
        // mais plus precis (au debut des iterations)
        // = expansion de la fonction hs_step0_F32(F0[level], F1r[level], h, w, alpha2, Fx[level], Fy[level], Ft[level], dU[level], dV[level]);

        iter = 0;

        hs_step0_F32 (F0[level], F1r[level], h, w, alpha2, Fx[level], Fy[level], Ft[level], dU[level], dV[level]);
        /*hs_calc_gradientXYT_F32(F0[level], F1r[level], h, w, Fx[level], Fy[level], Ft[level]);
        hs_calc_UV0_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]);/**/

        //puts("calcul avec propagation de (U,V) a chaque niveau et initialisation uniquement au plus bas niveau");

        /*if(level == nb_level - 1) {
            printf("level = %d == nb_level-1 = %d\n", level, nb_level-1);
            // si plus bas niveau de la pyramide, les vecteurs vitesse sont nuls: (U,V)=0 => (Um,Vm)=0
            // calcul des gradients avec comme seconde image, l'image initiale F1
            //hs_calc_gradientXYT_F32(F0[level], F1[level], h, w, Fx[level], Fy[level], Ft[level]); // grad(F0,F1) moins precis que grad(F0,F1r) = incoherent
            hs_calc_gradientXYT_F32(F0[level], F1r[level], h, w, Fx[level], Fy[level], Ft[level]);
            hs_calc_UV0_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]);

        } else {
            printf("level = %d\n", level);
            // sinon, les vecteurs vitesse ne sont pas nuls et sont dans (U,V)

            // calcul des gradients avec comme seconde image, l'image recalee F1r
            hs_calc_gradientXYT_F32(F0[level], F1r[level], h, w, Fx[level], Fy[level], Ft[level]);
            //hs_calc_UV0_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]);

            // faux mais converge mieux a 10 iterations (moins bonne convergence pour 5 iter)
            //hs_calc_UmVm_F32(U[level], V[level], h, w, Um[level], Vm[level]);
            //hs_calc_UV_F32(Fx[level], Fy[level], Ft[level], Um[level], Vm[level], h, w, alpha2, dU[level], dV[level]); // init de (dU,dV)

        }*/

        // -----------------
        // -- verbose MSE --
        // -----------------

        hs_add_UV_F32(U[level], V[level], dU[level], dV[level], h, w, Ut[level], Vt[level]); // addition dans Utotal Vtotal (temporaire)
        extendDataD2_F32(F1[level], h, w);
        interpCubic_F32(F1[level], Ut[level], Vt[level], h, w, F1t[level]);
        mse = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);
        printf("level = %d iter = %3d mse = %8.3f\n", level, iter, mse);
        //convert_f32matrix_ui8matrix_sat(F1t[level], 0, h-1, 0, w-1, I1r); // saturation necessaire
        //generate_path_filename_k_ndigit_l_extension(dst_path, "Ir_", level, nb_digit, iter, "pgm", complete_filename);
        //SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

        // ======================
        // == boucle iteration ==
        // ======================

        for(iter = 1; iter <= (nb_iter); iter++) {

            //printf("level = %d iter = %d\n", level, iter);

            // --------------------------------------------
            // -- 1step: calcul de (Um,Vm) et de (du,dv) --
            // --------------------------------------------

            // en deux fonctions
            hs_calc_UmVm_F32(dU[level], dV[level], h, w, Um[level], Vm[level]); // OK
            hs_calc_UV_F32(Fx[level], Fy[level], Ft[level], Um[level], Vm[level], h, w, alpha2, dU[level], dV[level]); //OK

            // en une fonction en appelant deux
            //hs_1step_F32(dU[level], dV[level], Fx[level], Fy[level], Ft[level], h, w, alpha2, Um[level], Vm[level], dU[level], dV[level]);

            // en une fonction calculant directement dU,dV sans passer par Um,Vm
            //hs_calc_UV_direct_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]); //OK 2017
            //hs_calc_UV_fusion_fast_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]); // version line pour pipeline futur
            //hs_calc_UV_red_fusion_fast_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]); // version line pour pipeline futur


            // -----------------
            // -- verbose MSE --
            // -----------------
            hs_add_UV_F32(U[level], V[level], dU[level], dV[level], h, w, Ut[level], Vt[level]); // (Ut,Vt) = (U,V) + (du,dv)
            extendDataD2_F32(F1[level], h, w);
            interpCubic_F32(F1[level], Ut[level], Vt[level], h, w, F1t[level]);
            mse = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);
            printf("level = %d iter = %3d mse = %8.3f\n", level, iter, mse);

            // ------------------
            // -- save to disk --
            // ------------------
            //convert_f32matrix_ui8matrix_sat(F1t[level], 0, h-1, 0, w-1, I1r); // saturation necessaire
            //generate_path_filename_k_ndigit_l_extension(dst_path, "Ir_", level, nb_digit, iter, "pgm", complete_filename);
            //SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

            if(level == nb_level - 1) {
                // generate_path_filename_k_ndigit_l_extension(dst_path, "dU_", level, nb_digit, iter, "txt", complete_filename); write_f32matrix(dU[level], 0, h-1, 0, w-1, "%8.3f", complete_filename);
                // generate_path_filename_k_ndigit_l_extension(dst_path, "dV_", level, nb_digit, iter, "txt", complete_filename); write_f32matrix(dV[level], 0, h-1, 0, w-1, "%8.3f", complete_filename);*/
            }

        } // iter

        // --------------
        // -- epilogue --
        // --------------

        hs_accumulate_UV_F32(dU[level], dV[level], h, w, U[level], V[level]); // (U,V) += (dU,dV)

        // -----------------
        // -- verbose MSE --
        // -----------------
        interpCubic_F32(F1[level], U[level], V[level], h, w, F1t[level]);
        mse = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);
        printf("Level = %d Iter = %3d MSE = %8.3f\n", level, nb_iter+1, mse);
        //convert_f32matrix_ui8matrix_sat(F1r[level], 0, h-1, 0, w-1, I1r); // saturation necessaire
        //generate_path_filename_k_ndigit_l_extension(dst_path, "I1r_", level, nb_digit, nb_iter+1, "pgm", complete_filename);
        //SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

        if(level) {
            resample_F32(U[level], h, w, U[level-1]);
            resample_F32(V[level], h, w, V[level-1]);

            // -- debug F16: OK
            /*if(level == nb_level - 1) {
             generate_path_filename_k_ndigit_extension(dst_path, "Uup_", level, nb_digit, "txt", complete_filename); write_f32matrix(U[level-1], 0, 2*h-1, 0, 2*w-1, "%8.3f", complete_filename);
             generate_path_filename_k_ndigit_extension(dst_path, "Vup_", level, nb_digit, "txt", complete_filename); write_f32matrix(V[level-1], 0, 2*h-1, 0, 2*w-1, "%8.3f", complete_filename);
             }*/
        }

    } // level
end:

    // free

    free_f32matrix(K, -2, +2, -2, +2);

    free_ui8matrix(I0,  i0, i1, j0, j1);
    free_ui8matrix(I1,  i0, i1, j0, j1);
    free_ui8matrix(I1r, i0, i1, j0, j1);

    free_f32pyramid(F0,  nb_level, height, width, border);
    free_f32pyramid(F1,  nb_level, height, width, border);
    free_f32pyramid(F1r, nb_level, height, width, border);
    free_f32pyramid(F1t, nb_level, height, width, border);

    free_f32pyramid(Fx, nb_level, height, width, border);
    free_f32pyramid(Fy, nb_level, height, width, border);
    free_f32pyramid(Ft, nb_level, height, width, border);

    free_f32pyramid(dU, nb_level, height, width, border);
    free_f32pyramid(dV, nb_level, height, width, border);

    free_f32pyramid(Um, nb_level, height, width, border);
    free_f32pyramid(Vm, nb_level, height, width, border);

    free_f32pyramid(U, nb_level, height, width, border);
    free_f32pyramid(V, nb_level, height, width, border);

    free_f32pyramid(Ut, nb_level, height, width, 0);
    free_f32pyramid(Vt, nb_level, height, width, 0);

    free_f32matrix(W, 0, height-1, 0, width-1);

    return 0;
}
// -------------------------------------------------------
int HS_test_hierarchique1_auto_F32(int argc, char *argv[])
// -------------------------------------------------------
{
    // HS_test_hierarchique1_F32 avec auto-tuning

    int h, w, height, width;  // image size
    int i0, i1, j0, j1; // image intervals
    int nb_digit = 2;
    int iter,  nb_iter  = 200;
    int level, nb_level = 3;
    int border = 3; // for interpolation

    float32 alpha2 = 0.005f;

    double mse;

    char *src_path = "/nfs/home/sasl/eleves/main/3800549/Bureau/OpticalFlow_HS/OpticalFlow_HS/config/";
    char *dst_path = "/Users/lacas/Code/OpticalFlow/resultHierarchique/";

    char *filename0;
    char *filename1;

    char complete_filename [1024];
    char complete_filename0[1024];
    char complete_filename1[1024];

    float32 *Alpha2;
    uint8   **I0, **I1, **I1r;
    float32 **K; // Burt Kernel

    float32 ***F0, ***F1, ***F1r, ***F1t;         // pair of images (t), (t+1), reconstruction, reconstruction temporaire
    float32 ***Fx, ***Fy, ***Ft;                  // {x, y, t} gradient
    float32 ***dU, ***dV;
    float32 ***Um, ***Vm;
    float32 ***U, ***V;
    float32 ***Ut, ***Vt;
    float32 **W;

    puts("---------------------------------");
    puts("--- HS_test_hierarchique1_F32 ---");
    puts("---------------------------------");

    // pour la paire d'images coastguard060 + coastguard061, 3 niveaux, 100 iters
    //alpha2 = 0.100f; // mse(100)=52.81
    //alpha2 = 0.050f; // mse(100)=43.54
    //alpha2 = 0.020f; // mse(100)=36.42
    //alpha2 = 0.010f; // mse(100)=33.95
    //alpha2 = 0.008f; // mse(100)=33.58
    //alpha2 = 0.007f; // mse(100)=33.46
    //alpha2 = 0.006f; // mse(100)=33.43
    //alpha2 = 0.005f; // mse(100)=33.53
    //alpha2 = 0.004f; // mse(100)=33.86
    //alpha2 = 0.003f; // mse(100)=35.42
    //alpha2 = 0.002f; // mse(100)=38.24
    //alpha2 = 0.001f; // mse(100)=43.77

    // pour la paire d'images coastguard150 + coastguard151, 3 niveaux, 100 iters

    //alpha2 = 0.100f; // mse = 46.64
    //alpha2 = 0.050f; // mse = 35.87
    //alpha2 = 0.020f; // mse = 25.14
    //alpha2 = 0.010f; // mse = 20.31
    //alpha2 = 0.008f; // mse = 19.34
    //alpha2 = 0.007f; // mse = 18.87
    //alpha2 = 0.006f; // mse = 18.43 <- best de (060,061)
    //alpha2 = 0.005f; // mse = 18.06
    //alpha2 = 0.004f; // mse = 17.81
    //alpha2 = 0.003f; // mse = 17.80 <- best
    //alpha2 = 0.002f; // mse = 18.31
    //alpha2 = 0.001f; // mse = 20.47

    alpha2 = 0.006f;
    alpha2 = 0.003f;

    nb_level = 3;
    nb_iter  = 20;

    //filename0 = "I0.pgm"; filename1 = "I1.pgm"; //synthetic cosine/sine
    //filename0 = "c60.pgm"; filename1 = "c60t.pgm";
    //filename0 = "meteor1_0495.pgm"; filename1 = "meteor1_0496.pgm";

    //filename0 = "coastguard060.pgm"; filename1 = "coastguard060t.pgm";
    //filename0 = "coastguard060.pgm"; filename1 = "coastguard061.pgm";

    //filename0 = "coastguard060.pgm"; filename1 = "coastguard060t.pgm";
    //filename0 = "coastguard060.pgm"; filename1 = "coastguard060q.pgm";

    //filename0 = "coastguard149.pgm"; filename1 = "coastguard150.pgm";
    filename0 = "coastguard150.pgm"; filename1 = "coastguard151.pgm"; // descend plus bas

    generate_path_filename(src_path, filename0, complete_filename0);
    generate_path_filename(src_path, filename1, complete_filename1);

    // load & alloc
    I0 = LoadPGM_ui8matrix(complete_filename0, &i0, &i1, &j0, &j1);
    I1 = LoadPGM_ui8matrix(complete_filename1, &i0, &i1, &j0, &j1);

    // full size
    height = i1 - i0 + 1;
    width  = j1 - j0 + 1;

    printf("size = %d x %d\n", width, height);
    printf("%s vs %s\n", filename0, filename1);
    printf("nb_level = %d\n", nb_level);
    printf("nb_iter = %d\n", nb_iter);

    // -----------
    // -- alloc --
    // -----------

    Alpha2 = f32vector(0, nb_level);
    I1r  = ui8matrix(0, height-1, 0, width-1);

    K = f32matrix(-2, +2, -2, +2); // 5x5 kernel for Burt Pyramid

    F0  = f32pyramid(nb_level, height, width, border);
    F1  = f32pyramid(nb_level, height, width, border);
    F1r = f32pyramid(nb_level, height, width, border);
    F1t = f32pyramid(nb_level, height, width, border);

    Fx = f32pyramid(nb_level, height, width, border);
    Fy = f32pyramid(nb_level, height, width, border);
    Ft = f32pyramid(nb_level, height, width, border);

    dU = f32pyramid(nb_level, height, width, border);
    dV = f32pyramid(nb_level, height, width, border);

    Um = f32pyramid(nb_level, height, width, border);
    Vm = f32pyramid(nb_level, height, width, border);

    U = f32pyramid(nb_level, height, width, border);
    V = f32pyramid(nb_level, height, width, border);

    Ut = f32pyramid(nb_level, height, width, 0);
    Vt = f32pyramid(nb_level, height, width, 0);

    W = f32matrix(0, height-1, 0, width-1);

    // ----------
    // -- init --
    // ----------

    zero_f32pyramid(F0,  nb_level, height, width, border);
    zero_f32pyramid(F1,  nb_level, height, width, border);
    zero_f32pyramid(F1r, nb_level, height, width, border);
    zero_f32pyramid(F1t, nb_level, height, width, border);

    zero_f32pyramid(Fx, nb_level, height, width, border);
    zero_f32pyramid(Fy, nb_level, height, width, border);
    zero_f32pyramid(Ft, nb_level, height, width, border);

    zero_f32pyramid(dU, nb_level, height, width, border);
    zero_f32pyramid(dV, nb_level, height, width, border);

    zero_f32pyramid(Um, nb_level, height, width, border);
    zero_f32pyramid(Vm, nb_level, height, width, border);

    zero_f32pyramid(U, nb_level, height, width, border);
    zero_f32pyramid(V, nb_level, height, width, border);

    zero_f32pyramid(Ut, nb_level, height, width, 0);
    zero_f32pyramid(Vt, nb_level, height, width, 0);

    zero_f32matrix(W, 0, height-1, 0, width-1);

    // convert avec division par 255: [0:255] -> [0:1]
    //convert_ui8matrix_f32matrix_quantif(I0, i0, i1, j0, j1, 255, E0);  // attention: pas 256
    //convert_ui8matrix_f32matrix_quantif(I1, i0, i1, j0, j1, 255, E1);

    // -------------------------------- //
    // -- construction des pyramides -- //
    // -------------------------------- //

    for(level = 0; level < nb_level; level++) { Alpha2[level] = alpha2; }

    alpha2 = 0.0020f;
    alpha2 = 0.0019f;
    Alpha2[2] = alpha2 * 255.f * 255.f;

    alpha2 = 0.0019f;
    alpha2 = 0.0018f;
    alpha2 = 0.0017f;
    alpha2 = 0.0016f;
    alpha2 = 0.0015f;
    alpha2 = 0.0014f;
    alpha2 = 0.0013f;
    alpha2 = 0.0012f;
    alpha2 = 0.0011f;
    alpha2 = 0.0010f;
    //alpha2 = 0.0009f;
    //alpha2 = 0.0008f;
    //alpha2 = 0.0005f;
    Alpha2[1] = alpha2 * 255.f * 255.f;

    alpha2 = 0.0020f;
    alpha2 = 0.0030f;
    alpha2 = 0.0025f;
    alpha2 = 0.0024f;
    //alpha2 = 0.0005f;
    Alpha2[0] = alpha2 * 255.f * 255.f;

    initBurtKernel_F32(0.4f, K);

    // prologue d'init

    level = 0;
    h = height >> level;
    w = width  >> level;

    convert_ui8matrix_f32matrix(I0, 0, h-1, 0, w-1, F0[level]);
    convert_ui8matrix_f32matrix(I1, 0, h-1, 0, w-1, F1[level]);

    // pour debug
    //convert_f32matrix_ui8matrix(F0[level], 0, h-1, 0, w-1, I0);
    //convert_f32matrix_ui8matrix(F1[level], 0, h-1, 0, w-1, I1);
    //generate_path_filename_k_ndigit_extension(dst_path, "I0_", level, nb_digit, "pgm", complete_filename0);
    //generate_path_filename_k_ndigit_extension(dst_path, "I1_", level, nb_digit, "pgm", complete_filename1);
    //SavePGM_ui8matrix(I0, 0, h-1, 0, w-1, complete_filename0);
    //SavePGM_ui8matrix(I1, 0, h-1, 0, w-1, complete_filename1);

    // -------------------
    // -- boucle d'init --
    // -------------------

    for(level = 1; level <= nb_level - 1 ; level++) {

        h = height >> level;
        w = width  >> level;

        // Burt pyramid
        convolution5x5_decimation_F32_F32_F32(F0[level-1], 2*h, 2*w, K, F0[level]);
        convolution5x5_decimation_F32_F32_F32(F1[level-1], 2*h, 2*w, K, F1[level]);

        // pour debug
        /*convert_f32matrix_ui8matrix(F0[level], 0, h-1, 0, w-1, I0);
         convert_f32matrix_ui8matrix(F1[level], 0, h-1, 0, w-1, I1);
         generate_path_filename_k_ndigit_extension(dst_path, "I0_", level, nb_digit, "pgm", complete_filename0);
         generate_path_filename_k_ndigit_extension(dst_path, "I1_", level, nb_digit, "pgm", complete_filename1);
         SavePGM_ui8matrix(I0, 0, h-1, 0, w-1, complete_filename0);
         SavePGM_ui8matrix(I1, 0, h-1, 0, w-1, complete_filename1);/**/

    } // level

    // ---------------- //
    // -- traitement -- //
    // ---------------- //

    // erreur initiale (full res)
    mse = of_calcMeanSquareError_F32(F0[0], F1[0], h, w);
    printf("level = %d iter = %3d mse = %8.4f\n", -1, -1, mse);

    // =====================================
    // == boucle pyramidale sur le niveau ==
    // =====================================

    for(level = nb_level - 1; level >= 0; level--) {

        iter = -1;
        printf("-----------------------------------\n");
        printf("level = %d", level);
        printf("  Alpha2[%d] = %.4f\n", level, Alpha2[level]/(255*255));
        h = height >> level;
        w = width  >> level;
        i0 = 0; i1 = h-1;
        j0 = 0; j1 = w-1;

        // erreur initiale au niveau courant, avant le recallage
        mse = of_calcMeanSquareError_F32(F0[level], F1[level], h, w);
        printf("level = %d iter = %3d mse = %8.4f avant warp initial\n", level, iter, mse);

        // inutile normalement car write before read
        //zero_f32matrix(dU[level], 0-border, h-1+border, 0-border, w-1+border);
        //zero_f32matrix(dV[level], 0-border, h-1+border, 0-border, w-1+border);

        // ----------------------------------
        // -- recalage de debut de niveau --
        // ----------------------------------

        // recalage / compensation du mouvement avec les vecteurs vitesses calcules et upsamples du niveau precedent

        extendDataD2_F32(F1[level], h, w);
        //if(level < nb_level-1) { interpCubic_F32(F1[level], U[level], V[level], h, w, F1r[level]); }

        interpCubic_F32(F1[level], U[level], V[level], h, w, F1r[level]);

        mse = of_calcMeanSquareError_F32(F0[level], F1r[level], h, w);
        printf("level = %d iter = %3d mse = %8.4f apres warp initial\n", level, iter, mse);

        // iteration initiale avec (U,V)=0 => (Um,Vm)=0
        /*if(level == nb_level-1) {
         hs_step0_F32(F0[level], F1[level], h, w, alpha2, Fx[level], Fy[level], Ft[level], dU[level], dV[level]);
         }*/

        // ----------------------------------------------------------------------------------
        // -- step0: calcul des gradients (Ix,Iy,It) et init des vecteurs vitesses (dU,dV) --
        // ----------------------------------------------------------------------------------

        // naivement (faux)
        //put("calcul avec reset de (U,V) a chaque niveau et initialisation via derivees");
        // mais plus precis (au debut des iterations)
        // = expansion de la fonction hs_step0_F32(F0[level], F1r[level], h, w, alpha2, Fx[level], Fy[level], Ft[level], dU[level], dV[level]);

        iter = 0;

        //hs_step0_F32 (F0[level], F1r[level], h, w, alpha2, Fx[level], Fy[level], Ft[level], dU[level], dV[level]);
        hs_step0_F32 (F0[level], F1r[level], h, w, Alpha2[level], Fx[level], Fy[level], Ft[level], dU[level], dV[level]);
        /*hs_calc_gradientXYT_F32(F0[level], F1r[level], h, w, Fx[level], Fy[level], Ft[level]);
         hs_calc_UV0_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]);/**/

        //puts("calcul avec propagation de (U,V) a chaque niveau et initialisation uniquement au plus bas niveau");

        /*if(level == nb_level - 1) {
         printf("level = %d == nb_level-1 = %d\n", level, nb_level-1);
         // si plus bas niveau de la pyramide, les vecteurs vitesse sont nuls: (U,V)=0 => (Um,Vm)=0
         // calcul des gradients avec comme seconde image, l'image initiale F1
         //hs_calc_gradientXYT_F32(F0[level], F1[level], h, w, Fx[level], Fy[level], Ft[level]); // grad(F0,F1) moins precis que grad(F0,F1r) = incoherent
         hs_calc_gradientXYT_F32(F0[level], F1r[level], h, w, Fx[level], Fy[level], Ft[level]);
         hs_calc_UV0_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]);

         } else {
         printf("level = %d\n", level);
         // sinon, les vecteurs vitesse ne sont pas nuls et sont dans (U,V)

         // calcul des gradients avec comme seconde image, l'image recalee F1r
         hs_calc_gradientXYT_F32(F0[level], F1r[level], h, w, Fx[level], Fy[level], Ft[level]);
         //hs_calc_UV0_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]);

         // faux mais converge mieux a 10 iterations (moins bonne convergence pour 5 iter)
         //hs_calc_UmVm_F32(U[level], V[level], h, w, Um[level], Vm[level]);
         //hs_calc_UV_F32(Fx[level], Fy[level], Ft[level], Um[level], Vm[level], h, w, alpha2, dU[level], dV[level]); // init de (dU,dV)

         }*/

        // -----------------
        // -- verbose MSE --
        // -----------------

        hs_add_UV_F32(U[level], V[level], dU[level], dV[level], h, w, Ut[level], Vt[level]); // addition dans Utotal Vtotal (temporaire)
        extendDataD2_F32(F1[level], h, w);
        interpCubic_F32(F1[level], Ut[level], Vt[level], h, w, F1t[level]);
        mse = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);
        printf("level = %d iter = %3d mse = %8.4f\n", level, iter, mse);
        //convert_f32matrix_ui8matrix_sat(F1t[level], 0, h-1, 0, w-1, I1r); // saturation necessaire
        //generate_path_filename_k_ndigit_l_extension(dst_path, "Ir_", level, nb_digit, iter, "pgm", complete_filename);
        //SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

        // ======================
        // == boucle iteration ==
        // ======================

        for(iter = 1; iter <= (nb_iter); iter++) {

            //printf("level = %d iter = %d\n", level, iter);

            // --------------------------------------------
            // -- 1step: calcul de (Um,Vm) et de (du,dv) --
            // --------------------------------------------

            // en deux fonctions
            hs_calc_UmVm_F32(dU[level], dV[level], h, w, Um[level], Vm[level]); // OK
            //hs_calc_UV_F32(Fx[level], Fy[level], Ft[level], Um[level], Vm[level], h, w, alpha2, dU[level], dV[level]); //OK
            hs_calc_UV_F32(Fx[level], Fy[level], Ft[level], Um[level], Vm[level], h, w, Alpha2[level], dU[level], dV[level]); //OK

            // en une fonction en appelant deux
            //hs_1step_F32(dU[level], dV[level], Fx[level], Fy[level], Ft[level], h, w, alpha2, Um[level], Vm[level], dU[level], dV[level]);

            // en une fonction calculant directement dU,dV sans passer par Um,Vm
            //hs_calc_UV_direct_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]); //OK 2017
            //hs_calc_UV_fusion_fast_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]); // version line pour pipeline futur
            //hs_calc_UV_red_fusion_fast_F32(Fx[level], Fy[level], Ft[level], h, w, alpha2, dU[level], dV[level]); // version line pour pipeline futur


            // -----------------
            // -- verbose MSE --
            // -----------------
            hs_add_UV_F32(U[level], V[level], dU[level], dV[level], h, w, Ut[level], Vt[level]); // (Ut,Vt) = (U,V) + (du,dv)
            extendDataD2_F32(F1[level], h, w);
            interpCubic_F32(F1[level], Ut[level], Vt[level], h, w, F1t[level]);
            mse = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);
            printf("level = %d iter = %3d mse = %8.4f\n", level, iter, mse);

        } // iter

        // --------------
        // -- epilogue --
        // --------------

        hs_accumulate_UV_F32(dU[level], dV[level], h, w, U[level], V[level]); // (U,V) += (dU,dV)

        // -----------------
        // -- verbose MSE --
        // -----------------
        interpCubic_F32(F1[level], U[level], V[level], h, w, F1t[level]);
        mse = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);
        printf("Level = %d Iter = %3d MSE = %8.4f\n", level, nb_iter+1, mse);

        if(level) {
            resample_F32(U[level], h, w, U[level-1]);
            resample_F32(V[level], h, w, V[level-1]);

        }

    } // level
end:

    // free

    free_f32vector(Alpha2, 0, nb_level);

    free_ui8matrix(I0,  i0, i1, j0, j1);
    free_ui8matrix(I1,  i0, i1, j0, j1);
    free_ui8matrix(I1r, i0, i1, j0, j1);

    free_f32pyramid(F0,  nb_level, height, width, border);
    free_f32pyramid(F1,  nb_level, height, width, border);
    free_f32pyramid(F1r, nb_level, height, width, border);
    free_f32pyramid(F1t, nb_level, height, width, border);

    free_f32pyramid(Fx, nb_level, height, width, border);
    free_f32pyramid(Fy, nb_level, height, width, border);
    free_f32pyramid(Ft, nb_level, height, width, border);

    free_f32pyramid(dU, nb_level, height, width, border);
    free_f32pyramid(dV, nb_level, height, width, border);

    free_f32pyramid(Um, nb_level, height, width, border);
    free_f32pyramid(Vm, nb_level, height, width, border);

    free_f32pyramid(U, nb_level, height, width, border);
    free_f32pyramid(V, nb_level, height, width, border);

    free_f32pyramid(Ut, nb_level, height, width, 0);
    free_f32pyramid(Vt, nb_level, height, width, 0);

    free_f32matrix(W, 0, height-1, 0, width-1);

    return 0;
}
// -------------------------------------------------
int HS_test_hierarchique1_FQ(int argc, char *argv[])
// -------------------------------------------------
{
#ifdef ENABLE_FIXED_POINT
    int h, w, height, width;  // image size
    int i0, i1, j0, j1; // image intervals
    int nb_digit = 2;
    int iter, nb_iter;
    int level, nb_level;
    int border = 3; // for interpolation

    int q_image, q_derivative, q_speed, q_interp, q_kernel; // speed quantification
    int Q_image, Q_derivative, Q_speed, Q_interp, Q_kernel;

    float32 alpha2; // valeur F32 pour des images de dynamique 1: [0..1] (et pas sur 0..255)
    float32 alpha2_f; // pour des images F32 quantifiees
    sint32  alpha2_q; // pour des images I16 quantifiees

    double mse_f, mse_q, mse_m, mse_c;

    char *src_path = "/nfs/home/sasl/eleves/main/3800549/Bureau/OpticalFlow_HS/OpticalFlow_HS/config/";
    char *dst_path = "/Users/lacas/Code/OpticalFlow/resultHierarchique/";

    char *filename0;
    char *filename1;

    char complete_filename [1024];
    char complete_filename0[1024];
    char complete_filename1[1024];

    uint8   **I0, **I1, **I1r;

    // Burt Kernel
    float32 **K_f;
    uint16 **K_q;

    float32 ***F0, ***F1, ***F1r, ***F1t;
    float32 ***Fx, ***Fy, ***Ft;
    float32 ***dU_f, ***dV_f;
    float32 ***Um_f, ***Vm_f;
    float32 ***U_f, ***V_f;
    float32 ***Ut_f, ***Vt_f;

    uint16 ***E0, ***E1;

    //uint16 ***E1t; // temporaire

    // Qstd
    uint16 ***E1r_q, ***E1t_q;
    sint16 ***Ex_q, ***Ey_q, ***Et_q;
    sint32 ***dU_q, ***dV_q;
    sint32 ***Um_q, ***Vm_q;
    sint32 ***U_q, ***V_q;
    sint32 ***Ut_q, ***Vt_q;

    // Qmax
    uint16 ***E1r_m, ***E1t_m;
    sint16 ***Ex_m, ***Ey_m, ***Et_m;
    sint32 ***dU_m, ***dV_m;
    sint32 ***Um_m, ***Vm_m;
    sint32 ***U_m, ***V_m;
    sint32 ***Ut_m, ***Vt_m;

    // conversion Q32 -> F32
    float32 ***dU_c, ***dV_c;
    float32  ***U_c,  ***V_c;
    float32  ***Ut_c, ***Vt_c;

    float32 **W;

    puts("---------------------------------");
    puts("--- HS_test_hierarchique1_F+Q ---");
    puts("---------------------------------");

    alpha2 = 0.635f; // mauvaise convergence
    alpha2 = 0.05f;  // bonne convergence
    alpha2 = 0.04f;
    alpha2 = 0.02f;
    alpha2 = 0.01f;
    //alpha2 = 0.005f;
    //alpha2 = 0.002f;
    //alpha2 = 0.001f;

    alpha2 = 0.006f;
    //alpha2 = 0.003f;

    nb_level = 3;
    nb_iter  = 10;

    nb_level = 3;
    nb_iter  = 20;

    // quantification 'basique'
    q_image      = 8;
    q_derivative = 8;
    q_speed      = 8;
    q_interp     = 8;
    q_kernel     = 8;

    // quantification 'embedded'
    /*q_image      = 8;
     q_derivative = 7;
     q_speed      = 4;
     q_interp     = 5;*/

    Q_image      = CALC_Q(q_image);
    Q_derivative = CALC_Q(q_derivative);
    Q_speed      = CALC_Q(q_speed);
    Q_interp     = CALC_Q(q_interp);
    Q_kernel     = CALC_Q(q_kernel);

    alpha2_f = alpha2 * Q_image * Q_image;   // en F32, images et derivees sont quantifiees sur le meme nombre de bits
    alpha2_q = (int) roundf(alpha2 * Q_derivative * Q_derivative); // en Q, images et derivves sont quantifies sur un nombre different de bits

    //filename0 = "I0.pgm"; filename1 = "I1.pgm"; //synthetic cosine/sine
    //filename0 = "c60.pgm"; filename1 = "c60t.pgm";
    //filename0 = "meteor1_0495.pgm"; filename1 = "meteor1_0496.pgm";

    //filename0 = "coastguard060.pgm"; filename1 = "coastguard060t.pgm";
    //filename0 = "coastguard060.pgm"; filename1 = "coastguard061.pgm";

    //filename0 = "coastguard060.pgm"; filename1 = "coastguard060t.pgm";
    filename0 = "coastguard060.pgm"; filename1 = "coastguard060q.pgm";

    //filename0 = "coastguard149.pgm"; filename1 = "coastguard150.pgm";
    filename0 = "coastguard150.pgm"; filename1 = "coastguard151.pgm"; // descend plus bas

    generate_path_filename(src_path, filename0, complete_filename0);
    generate_path_filename(src_path, filename1, complete_filename1);

    // load & alloc
    I0 = LoadPGM_ui8matrix(complete_filename0, &i0, &i1, &j0, &j1);
    I1 = LoadPGM_ui8matrix(complete_filename1, &i0, &i1, &j0, &j1);

    // full size
    height = i1 - i0 + 1;
    width  = j1 - j0 + 1;

    int i0d, i1d, j0d, j1d;

    i0d = 1*height/4; i1d = 3*height/4; j0d = 1*width/4; j1d = 1*width/4; // 1/4 .. 3/4
    i0d = 1*height/8; i1d = 2*height/8; j0d = 1*width/8; j1d = 2*width/8; // 1/8 .. 2/8
    i0d = 10; i1d = i0d + 7; j0d = 10; j1d = j0d + 7; //

    printf("size = %d x %d\n", width, height);
    printf("%s vs %s\n", filename0, filename1);
    printf("nb_level = %d\n", nb_level);
    printf("nb_iter = %d\n", nb_iter);

    // -----------
    // -- alloc --
    // -----------

    I1r = ui8matrix(0, height-1, 0, width-1);

    K_f = f32matrix( -2, +2, -2, +2); // 5x5 kernel for Burt Pyramid
    K_q = ui16matrix(-2, +2, -2, +2);

    F0  = f32pyramid(nb_level, height, width, border);
    F1  = f32pyramid(nb_level, height, width, border);
    F1r = f32pyramid(nb_level, height, width, border);
    F1t = f32pyramid(nb_level, height, width, border);

    Fx = f32pyramid(nb_level, height, width, border);
    Fy = f32pyramid(nb_level, height, width, border);
    Ft = f32pyramid(nb_level, height, width, border);

    dU_f = f32pyramid(nb_level, height, width, border);
    dV_f = f32pyramid(nb_level, height, width, border);

    Um_f = f32pyramid(nb_level, height, width, border);
    Vm_f = f32pyramid(nb_level, height, width, border);

    U_f = f32pyramid(nb_level, height, width, border);
    V_f = f32pyramid(nb_level, height, width, border);

    Ut_f = f32pyramid(nb_level, height, width, 0);
    Vt_f = f32pyramid(nb_level, height, width, 0);

    E0  = ui16pyramid(nb_level, height, width, border);
    E1  = ui16pyramid(nb_level, height, width, border);

    // Qstd
    E1r_q = ui16pyramid(nb_level, height, width, border);
    E1t_q = ui16pyramid(nb_level, height, width, border);

    Ex_q = si16pyramid(nb_level, height, width, border);
    Ey_q = si16pyramid(nb_level, height, width, border);
    Et_q = si16pyramid(nb_level, height, width, border);

    dU_q = si32pyramid(nb_level, height, width, border);
    dV_q = si32pyramid(nb_level, height, width, border);

    Um_q = si32pyramid(nb_level, height, width, border);
    Vm_q = si32pyramid(nb_level, height, width, border);

    U_q = si32pyramid(nb_level, height, width, border);
    V_q = si32pyramid(nb_level, height, width, border);

    Ut_q = si32pyramid(nb_level, height, width, 0);
    Vt_q = si32pyramid(nb_level, height, width, 0);

    // Qmax
    E1r_m = ui16pyramid(nb_level, height, width, border);
    E1t_m = ui16pyramid(nb_level, height, width, border);

    Ex_m = si16pyramid(nb_level, height, width, border);
    Ey_m = si16pyramid(nb_level, height, width, border);
    Et_m = si16pyramid(nb_level, height, width, border);

    dU_m = si32pyramid(nb_level, height, width, border);
    dV_m = si32pyramid(nb_level, height, width, border);

    Um_m = si32pyramid(nb_level, height, width, border);
    Vm_m = si32pyramid(nb_level, height, width, border);

    U_m = si32pyramid(nb_level, height, width, border);
    V_m = si32pyramid(nb_level, height, width, border);

    Ut_m = si32pyramid(nb_level, height, width, 0);
    Vt_m = si32pyramid(nb_level, height, width, 0);

    // convert
    dU_c = f32pyramid(nb_level, height, width, border);
    dV_c = f32pyramid(nb_level, height, width, border);
    U_c  = f32pyramid(nb_level, height, width, border);
    V_c  = f32pyramid(nb_level, height, width, border);
    Ut_c = f32pyramid(nb_level, height, width, 0);
    Vt_c = f32pyramid(nb_level, height, width, 0);

    W = f32matrix(0, height-1, 0, width-1);

    // ----------
    // -- init --
    // ----------

    // float
    zero_f32pyramid(F0,  nb_level, height, width, border);
    zero_f32pyramid(F1,  nb_level, height, width, border);
    zero_f32pyramid(F1r, nb_level, height, width, border);
    zero_f32pyramid(F1t, nb_level, height, width, border);

    zero_f32pyramid(Fx, nb_level, height, width, border);
    zero_f32pyramid(Fy, nb_level, height, width, border);
    zero_f32pyramid(Ft, nb_level, height, width, border);

    zero_f32pyramid(dU_f, nb_level, height, width, border);
    zero_f32pyramid(dV_f, nb_level, height, width, border);

    zero_f32pyramid(Um_f, nb_level, height, width, border);
    zero_f32pyramid(Vm_f, nb_level, height, width, border);

    zero_f32pyramid(U_f, nb_level, height, width, border);
    zero_f32pyramid(V_f, nb_level, height, width, border);

    zero_f32pyramid(Ut_f, nb_level, height, width, 0);
    zero_f32pyramid(Vt_f, nb_level, height, width, 0);

    // Qstd
    zero_ui16pyramid(E0,  nb_level, height, width, border);
    zero_ui16pyramid(E1,  nb_level, height, width, border);

    zero_ui16pyramid(E1r_q, nb_level, height, width, border);
    zero_ui16pyramid(E1t_q, nb_level, height, width, border);

    zero_si16pyramid(Ex_q, nb_level, height, width, border);
    zero_si16pyramid(Ey_q, nb_level, height, width, border);
    zero_si16pyramid(Et_q, nb_level, height, width, border);

    zero_si32pyramid(dU_q, nb_level, height, width, border);
    zero_si32pyramid(dV_q, nb_level, height, width, border);

    zero_si32pyramid(Um_q, nb_level, height, width, border);
    zero_si32pyramid(Vm_q, nb_level, height, width, border);

    zero_si32pyramid(U_q, nb_level, height, width, border);
    zero_si32pyramid(V_q, nb_level, height, width, border);

    zero_si32pyramid(Ut_q, nb_level, height, width, 0);
    zero_si32pyramid(Vt_q, nb_level, height, width, 0);

    // Qmax
    zero_ui16pyramid(E1r_m, nb_level, height, width, border);
    zero_ui16pyramid(E1t_m, nb_level, height, width, border);

    zero_si16pyramid(Ex_m, nb_level, height, width, border);
    zero_si16pyramid(Ey_m, nb_level, height, width, border);
    zero_si16pyramid(Et_m, nb_level, height, width, border);

    zero_si32pyramid(dU_m, nb_level, height, width, border);
    zero_si32pyramid(dV_m, nb_level, height, width, border);

    zero_si32pyramid(Um_m, nb_level, height, width, border);
    zero_si32pyramid(Vm_m, nb_level, height, width, border);

    zero_si32pyramid(U_m, nb_level, height, width, border);
    zero_si32pyramid(V_m, nb_level, height, width, border);

    zero_si32pyramid(Ut_m, nb_level, height, width, 0);
    zero_si32pyramid(Vt_m, nb_level, height, width, 0);

    // convert
    zero_f32pyramid(dU_c, nb_level, height, width, border);
    zero_f32pyramid(dV_c, nb_level, height, width, border);
    zero_f32pyramid(U_c, nb_level, height, width, border);
    zero_f32pyramid(V_c, nb_level, height, width, border);
    zero_f32pyramid(Ut_c, nb_level, height, width, 0);
    zero_f32pyramid(Vt_c, nb_level, height, width, 0);

    zero_f32matrix(W, 0, height-1, 0, width-1);

    // convert avec division par 255: [0:255] -> [0:1]
    //convert_ui8matrix_f32matrix_quantif(I0, i0, i1, j0, j1, 255, E0);  // attention: pas 256
    //convert_ui8matrix_f32matrix_quantif(I1, i0, i1, j0, j1, 255, E1);

    // -------------------------------- //
    // -- construction des pyramides -- //
    // -------------------------------- //

    initBurtKernel_F32(0.4f, K_f);  // F32
    quantifKernel5_q16(K_f, q_kernel, K_q);  // I16

    // prologue d'init

    level = 0;
    h = height >> level;
    w = width  >> level;

    convert_ui8matrix_f32matrix(I0, 0, h-1, 0, w-1, F0[level]); // F32
    convert_ui8matrix_f32matrix(I1, 0, h-1, 0, w-1, F1[level]);

    convert_ui8matrix_ui16matrix(I0, 0, h-1, 0, w-1, E0[level]); // I16
    convert_ui8matrix_ui16matrix(I1, 0, h-1, 0, w-1, E1[level]);

    // -------------------
    // -- boucle d'init --
    // -------------------

    for(level = 1; level <= nb_level - 1 ; level++) {

        h = height >> level;
        w = width  >> level;

        // Burt pyramid
        convolution5x5_decimation_F32_F32_F32(F0[level-1], 2*h, 2*w, K_f, F0[level]); // F32
        convolution5x5_decimation_F32_F32_F32(F1[level-1], 2*h, 2*w, K_f, F1[level]);

        convolution5x5_decimation_U16_I64_U16(E0[level-1], 2*h, 2*w, K_q, q_kernel, E0[level]); // I16
        convolution5x5_decimation_U16_I64_U16(E1[level-1], 2*h, 2*w, K_q, q_kernel, E1[level]);

        // pour debug
        /*convert_f32matrix_ui8matrix(F0[level], 0, h-1, 0, w-1, I0);
        convert_f32matrix_ui8matrix(F1[level], 0, h-1, 0, w-1, I1);
        generate_path_filename_k_ndigit_extension(dst_path, "F0_", level, nb_digit, "pgm", complete_filename0);
        generate_path_filename_k_ndigit_extension(dst_path, "F1_", level, nb_digit, "pgm", complete_filename1);
        SavePGM_ui8matrix(I0, 0, h-1, 0, w-1, complete_filename0);
        SavePGM_ui8matrix(I1, 0, h-1, 0, w-1, complete_filename1);

        convert_ui16matrix_ui8matrix(E0[level], 0, h-1, 0, w-1, I0);
        convert_ui16matrix_ui8matrix(E1[level], 0, h-1, 0, w-1, I1);
        generate_path_filename_k_ndigit_extension(dst_path, "E0_", level, nb_digit, "pgm", complete_filename0);
        generate_path_filename_k_ndigit_extension(dst_path, "E1_", level, nb_digit, "pgm", complete_filename1);
        SavePGM_ui8matrix(I0, 0, h-1, 0, w-1, complete_filename0);
        SavePGM_ui8matrix(I1, 0, h-1, 0, w-1, complete_filename1);*/

    } // level

    // ---------------- //
    // -- traitement -- //
    // ---------------- //

    // erreur initiale (full res)
    mse_f = of_calcMeanSquareError_F32(F0[0], F1[0], h, w);
    mse_q = of_calcMeanSquareError_U16(E0[0], E1[0], h, w);
    mse_m = of_calcMeanSquareError_U16(E0[0], E1[0], h, w);
    printf("level = %d iter = %3d  mse_f = %5.1f  mse_q = %5.1f  mse_m = %5.1f\n", -1, -1, mse_f, mse_q, mse_m);

    // =====================================
    // == boucle pyramidale sur le niveau ==
    // =====================================

    for(level = nb_level - 1; level >= 0; level--) {

        iter = -1;
        printf("-----------------------------------\n");
        printf("level = %d\n", level);
        h = height >> level;
        w = width  >> level;
        i0 = 0; i1 = h-1;
        j0 = 0; j1 = w-1;

        // erreur initiale au niveau courant, avant le recallage
        mse_f = of_calcMeanSquareError_F32(F0[level], F1[level], h, w);
        mse_q = of_calcMeanSquareError_U16(E0[level], E1[level], h, w);
        mse_m = of_calcMeanSquareError_U16(E0[level], E1[level], h, w);
        printf("level = %d iter = %3d  mse_f = %5.1f  mse_q = %5.1f  mse_m = %5.1f avant warp initial\n", level, iter, mse_f, mse_q, mse_m);

        // inutile normalement car write before read
        zero_f32matrix(dU_f[level], 0-border, h-1+border, 0-border, w-1+border);
        zero_f32matrix(dV_f[level], 0-border, h-1+border, 0-border, w-1+border);

        zero_si32matrix(dU_q[level], 0-border, h-1+border, 0-border, w-1+border);
        zero_si32matrix(dV_q[level], 0-border, h-1+border, 0-border, w-1+border);

        // ---------------------------------
        // -- recalage de debut de niveau --
        // ---------------------------------

        extendDataD2_F32(F1[level], h, w);
        extendDataD2_I16(E1[level], h, w);

        // conversion de debut de niveau (U_q,V_q) -> (U_c,V_c)
        convert_si32matrix_f32matrix_quantif(U_q[level], 0, h-1, 0, w-1, Q_speed, U_c[level]);
        convert_si32matrix_f32matrix_quantif(V_q[level], 0, h-1, 0, w-1, Q_speed, V_c[level]);

        interpCubic_F32    (F1[level],          U_f[level], V_f[level],          h, w, F1r[level]); mse_f = of_calcMeanSquareError_F32(F0[level], F1r[level], h, w);
        interpCubic_I16_I32(E1[level], q_image, U_q[level], V_q[level], q_speed, h, w, E1r_q[level]); mse_q = of_calcMeanSquareError_U16(E0[level], E1r_q[level], h, w);
        interpCubic_I16_I32(E1[level], q_image, U_m[level], V_m[level], q_speed, h, w, E1r_m[level]); mse_m = of_calcMeanSquareError_U16(E0[level], E1r_m[level], h, w);
        //interpCubic_I16_F32(E1[level], q_image, U_q[level], V_q[level], q_speed, h, w, E1r_q[level]); mse_q = of_calcMeanSquareError_U16(E0[level], E1r_q[level], h, w);
        interpCubic_F32    (F1[level],          U_c[level], V_c[level],          h, w, F1t[level]); mse_c = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);

        printf("level = %d iter = %3d  mse_f = %5.1f  mse_q = %5.1f  mse_m = %5.1f  mse_c = %5.1f  apres warp initial\n", level, 0, mse_f, mse_q, mse_m, mse_c);

        // ----------------------------------------------------------------------------------
        // -- step0: calcul des gradients (Ix,Iy,It) et init des vecteurs vitesses (dU,dV) --
        // ----------------------------------------------------------------------------------

        // iteration initiale avec (U,V)=0 => (Um,Vm)=0
        iter = 0;
        hs_step0_F32 (F0[level], F1r[level],           h, w, alpha2_f, Fx[level], Fy[level], Ft[level],               dU_f[level], dV_f[level]);
        hs_step0_Qstd(E0[level], E1r_q[level],  q_image, h, w, alpha2_q, Ex_q[level], Ey_q[level], Et_q[level], q_derivative, dU_q[level], dV_q[level], q_speed);
        hs_step0_Qmax(E0[level], E1r_m[level],  q_image, h, w, alpha2_q, Ex_m[level], Ey_m[level], Et_m[level], q_derivative, dU_m[level], dV_m[level], q_speed);
        //hs_step0_FQstd(E0[level], E1r[level],  q_image, h, w, alpha2_q, Ex[level], Ey[level], Et[level], q_derivative, dU_q[level], dV_q[level], q_speed);

        // conversion (dU_q,dV_q) -> (dU_c,dV_c)
        convert_si32matrix_f32matrix_quantif(dU_q[level], 0, h-1, 0, w-1, Q_speed, dU_c[level]);
        convert_si32matrix_f32matrix_quantif(dV_q[level], 0, h-1, 0, w-1, Q_speed, dV_c[level]);

        // -------------------
        // -- verbose MSE_F --
        // -------------------
        hs_add_UV_F32(U_f[level], V_f[level], dU_f[level], dV_f[level], h, w, Ut_f[level], Vt_f[level]); // addition dans Utotal Vtotal (temporaire)
        extendDataD2_F32(F1[level], h, w);
        interpCubic_F32(F1[level], Ut_f[level], Vt_f[level], h, w, F1t[level]);
        mse_f = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);

        // -------------------
        // -- verbose MSE_Q --
        // -------------------
        hs_add_UV_I32(U_q[level], V_q[level], dU_q[level], dV_q[level], h, w, Ut_q[level], Vt_q[level]); // addition dans Utotal Vtotal (temporaire)
        extendDataD2_I16(E1[level], h, w);
        interpCubic_I16_I32(E1[level], q_image, Ut_q[level], Vt_q[level], q_speed, h, w, E1t_q[level]);
        //interpCubic_I16_F32(E1[level], q_image, Ut_q[level], Vt_q[level], q_speed, h, w, E1t_q[level]);
        mse_q = of_calcMeanSquareError_U16(E0[level], E1t_q[level], h, w);

        // -------------------
        // -- verbose MSE_M --
        // -------------------
        hs_add_UV_I32(U_m[level], V_m[level], dU_m[level], dV_m[level], h, w, Ut_m[level], Vt_m[level]); // addition dans Utotal Vtotal (temporaire)
        extendDataD2_I16(E1[level], h, w);
        interpCubic_I16_I32(E1[level], q_image, Ut_q[level], Vt_q[level], q_speed, h, w, E1t_m[level]);
        //interpCubic_I16_F32(E1[level], q_image, Ut_m[level], Vt_m[level], q_speed, h, w, E1t_m[level]);
        mse_m = of_calcMeanSquareError_U16(E0[level], E1t_m[level], h, w);

        // -------------------
        // -- verbose MSE_C --
        // -------------------
        // version precision: addition dans F
        hs_add_UV_F32(U_c[level], V_c[level], dU_c[level], dV_c[level], h, w, Ut_c[level], Vt_c[level]); // addition dans Utotal Vtotal (temporaire)
        extendDataD2_F32(F1[level], h, w);
        interpCubic_F32(F1[level], Ut_c[level], Vt_c[level], h, w, F1t[level]);
        mse_c = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);

        //printf("level = %d iter = %3d  mse_f = %7.2f  mse_q = %7.2f apres warp initial\n", level, 0, mse_f, mse_q);
        printf("level = %d iter = %3d  mse_f = %5.1f  mse_q = %5.1f  mse_m = %5.1f  mse_c = %5.1f  apres warp initial\n", level, 0, mse_f, mse_q, mse_m, mse_c);

        // ----------------------
        // -- comparaison F32 I32
        // ----------------------
        //double mse_u = of_calcMeanSquareError_S32_F32(U_q[level], q_speed, U_f[level], h, w);
        //double mse_v = of_calcMeanSquareError_S32_F32(V_q[level], q_speed, V_f[level], h, w);

        //double mse_du = of_calcMeanSquareError_S32_F32(dU_q[level], q_speed, dU_f[level], h, w);
        //double mse_dv = of_calcMeanSquareError_S32_F32(dV_q[level], q_speed, dV_f[level], h, w);

        //printf("level = %d iter = %3d  mse_u = %10.6f  mse_v = %10.6f   mse_du = %10.6f  mse_db = %10.6f\n", level, 0, mse_u, mse_v, mse_du, mse_dv);

        //double mse_x = of_calcMeanSquareError_S16_F32(Ex[level], 0, Fx[level], h, w);
        //double mse_y = of_calcMeanSquareError_S16_F32(Ey[level], 0, Fy[level], h, w);
        //double mse_t = of_calcMeanSquareError_S16_F32(Et[level], 0, Ft[level], h, w);

        //printf("level = %d iter = %3d  mse_x = %10.6f  mse_y = %10.6f   mse_t = %10.6f\n", level, 0, mse_x, mse_y, mse_t);

        // ------------------
        // -- save to disk --
        // ------------------

        // debug
        //convert_f32matrix_ui8matrix_sat(F1r[level], 0, h-1, 0, w-1, I1r);
        //generate_path_filename_k_ndigit_l_extension(dst_path, "F1r_", level, nb_digit, iter, "pgm", complete_filename);
        //SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

        convert_ui16matrix_ui8matrix_sat(E1r_q[level], 0, h-1, 0, w-1, I1r);
        generate_path_filename_k_ndigit_l_extension(dst_path, "E1r_q_", level, nb_digit, iter, "pgm", complete_filename);
        SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

        convert_ui16matrix_ui8matrix_sat(E1r_m[level], 0, h-1, 0, w-1, I1r);
        generate_path_filename_k_ndigit_l_extension(dst_path, "E1r_m_", level, nb_digit, iter, "pgm", complete_filename);
        SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

        /*display_si16matrix_quantif(Ex[level], i0d, i1d, j0d, j1d, 1, "%8.3f", "Ex Q16");
        display_si16matrix_quantif(Ey[level], i0d, i1d, j0d, j1d, 1, "%8.3f", "Ey Q16");
        display_si16matrix_quantif(Et[level], i0d, i1d, j0d, j1d, 1, "%8.3f", "Et Q16");

        display_f32matrix(Fx[level], i0d, i1d, j0d, j1d, "%8.3f", "Fx F32");
        display_f32matrix(Fy[level], i0d, i1d, j0d, j1d, "%8.3f", "Fy F32");
        display_f32matrix(Ft[level], i0d, i1d, j0d, j1d, "%8.3f", "Ft F32");*/

        /*display_si32matrix_quantif(dU_q[level], i0d, i1d, j0d, j1d, Q_speed, "%8.3f", "du I32");
        display_si32matrix_quantif(dV_q[level], i0d, i1d, j0d, j1d, Q_speed, "%8.3f", "dv I32");

        display_f32matrix(dU_f[level], i0d, i1d, j0d, j1d, "%8.3f", "du F32");
        display_f32matrix(dV_f[level], i0d, i1d, j0d, j1d, "%8.3f", "dv F32");*/

        // ======================
        // == boucle iteration ==
        // ======================

        for(iter = 1; iter <= nb_iter; iter++) {

            //printf("level = %d iter = %d\n", level, iter);

            // --------------------------------------------
            // -- 1step: calcul de (Um,Vm) et de (du,dv) --
            // --------------------------------------------

            hs_calc_UmVm_F32(dU_f[level], dV_f[level], h, w, Um_f[level], Vm_f[level]);
            hs_calc_UV_F32(Fx[level], Fy[level], Ft[level], Um_f[level], Vm_f[level], h, w, alpha2_f, dU_f[level], dV_f[level]);

            hs_calc_UmVm_Qstd(dU_q[level], dV_q[level], q_speed, h, w, Um_q[level], Vm_q[level]);
            hs_calc_UV_Qstd(Ex_q[level], Ey_q[level], Et_q[level], q_derivative, Um_q[level], Vm_q[level], q_speed, h, w, alpha2_q, dU_q[level], dV_q[level]);

            hs_calc_UmVm_Qmax(dU_m[level], dV_m[level], h, w, Um_m[level], Vm_m[level]);
            hs_calc_UV_Qmax(Ex_m[level], Ey_m[level], Et_m[level], q_derivative, Um_m[level], Vm_m[level], q_speed, h, w, alpha2_q, dU_m[level], dV_m[level]);

            //hs_calc_UmVm_FQstd(dU_q[level], dV_q[level], q_speed, h, w, Um_q[level], Vm_q[level]);
            //hs_calc_UV_FQstd(Ex[level], Ey[level], Et[level], q_derivative, Um_q[level], Vm_q[level], q_speed, h, w, alpha2_q, dU_q[level], dV_q[level]);

            // -------------------
            // -- verbose MSE_F --
            // -------------------
            hs_add_UV_F32(U_f[level], V_f[level], dU_f[level], dV_f[level], h, w, Ut_f[level], Vt_f[level]); // addition dans Utotal Vtotal (temporaire)
            extendDataD2_F32(F1[level], h, w);
            interpCubic_F32(F1[level], Ut_f[level], Vt_f[level], h, w, F1t[level]);
            mse_f = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);

            // -------------------
            // -- verbose MSE_Q --
            // -------------------
            hs_add_UV_I32(U_q[level], V_q[level], dU_q[level], dV_q[level], h, w, Ut_q[level], Vt_q[level]); // addition dans Utotal Vtotal (temporaire)
            extendDataD2_I16(E1[level], h, w);
            interpCubic_I16_I32(E1[level], q_image, Ut_q[level], Vt_q[level], q_speed, h, w, E1t_q[level]);
            //interpCubic_I16_F32(E1[level], q_image, Ut_q[level], Vt_q[level], q_speed, h, w, E1t_q[level]);
            mse_q = of_calcMeanSquareError_U16(E0[level], E1t_q[level], h, w);

            // -------------------
            // -- verbose MSE_M --
            // -------------------
            hs_add_UV_I32(U_m[level], V_m[level], dU_m[level], dV_m[level], h, w, Ut_m[level], Vt_m[level]); // addition dans Utotal Vtotal (temporaire)
            extendDataD2_I16(E1[level], h, w);
            interpCubic_I16_I32(E1[level], q_image, Ut_m[level], Vt_m[level], q_speed, h, w, E1t_m[level]);
            //interpCubic_I16_F32(E1[level], q_image, Ut_q[level], Vt_q[level], q_speed, h, w, E1t_m[level]);
            mse_m = of_calcMeanSquareError_U16(E0[level], E1t_m[level], h, w);

            // -------------------
            // -- verbose MSE_C --
            // -------------------
            convert_si32matrix_f32matrix_quantif(dU_q[level], 0, h-1, 0, w-1, Q_speed, dU_c[level]); // conversion (dU_q,dV_q) -> (dU_c,dV_c)
            convert_si32matrix_f32matrix_quantif(dV_q[level], 0, h-1, 0, w-1, Q_speed, dV_c[level]);
            hs_add_UV_F32(U_c[level], V_c[level], dU_c[level], dV_c[level], h, w, Ut_c[level], Vt_c[level]); // addition dans Utotal Vtotal (temporaire)
            extendDataD2_F32(F1[level], h, w);
            interpCubic_F32(F1[level], Ut_c[level], Vt_c[level], h, w, F1t[level]);
            mse_c = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);

            // ------------------
            // -- save to disk --
            // ------------------
            convert_ui16matrix_ui8matrix_sat(E1t_q[level], 0, h-1, 0, w-1, I1r); // saturation necessaire
            generate_path_filename_k_ndigit_l_extension(dst_path, "E1r_q_", level, nb_digit, iter, "pgm", complete_filename);
            SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

            convert_ui16matrix_ui8matrix_sat(E1t_m[level], 0, h-1, 0, w-1, I1r); // saturation necessaire
            generate_path_filename_k_ndigit_l_extension(dst_path, "E1r_m_", level, nb_digit, iter, "pgm", complete_filename);
            SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

            //printf("level = %d iter = %3d  mse_f = %7.2f  mse_q = %7.2f\n", level, iter, mse_f, mse_q);
            printf("level = %d iter = %3d  mse_f = %5.1f  mse_q = %5.1f  mse_m = %5.1f  mse_c = %5.1f\n", level, iter, mse_f, mse_q, mse_m, mse_c);

            if(level == nb_level - 1) {
                //generate_path_filename_k_ndigit_l_extension(dst_path, "dU_", level, nb_digit, iter, "txt", complete_filename); write_f32matrix(dU[level], 0, h-1, 0, w-1, "%8.3f", complete_filename);
                //generate_path_filename_k_ndigit_l_extension(dst_path, "dV_", level, nb_digit, iter, "txt", complete_filename);  write_f32matrix(dV[level], 0, h-1, 0, w-1, "%8.3f", complete_filename);
            }

        } // iter

        // --------------
        // -- epilogue --
        // --------------

        hs_accumulate_UV_F32(dU_f[level], dV_f[level], h, w, U_f[level], V_f[level]);
        hs_accumulate_UV_I32(dU_q[level], dV_q[level], h, w, U_q[level], V_q[level]);

        // conversion (dU_q,dV_q) -> (dU_c,dV_c)
        convert_si32matrix_f32matrix_quantif(U_q[level], 0, h-1, 0, w-1, Q_speed, U_c[level]);
        convert_si32matrix_f32matrix_quantif(V_q[level], 0, h-1, 0, w-1, Q_speed, V_c[level]);

        interpCubic_F32    (F1[level],          U_f[level], V_f[level],          h, w, F1r[level]);
        interpCubic_I16_I32(E1[level], q_image, U_q[level], V_q[level], q_speed, h, w, E1r_q[level]);
        interpCubic_I16_I32(E1[level], q_image, U_q[level], V_q[level], q_speed, h, w, E1r_m[level]);
        interpCubic_F32    (F1[level],          U_c[level], V_c[level],          h, w, F1t[level]);

        mse_f = of_calcMeanSquareError_F32(F0[level], F1r[level], h, w);
        mse_q = of_calcMeanSquareError_U16(E0[level], E1r_q[level], h, w);
        mse_m = of_calcMeanSquareError_U16(E0[level], E1r_m[level], h, w);
        mse_c = of_calcMeanSquareError_F32(F0[level], F1t[level], h, w);

        //printf("level = %d iter = %3d  MSE_F = %5.1f  MSE_Q = %5.1f\n", level, nb_iter+1, mse_f, mse_q);
        printf("level = %d iter = %3d  MSE_F = %5.1f  MSE_Q = %5.1f  MSE_M = %5.1f  MSE_C = %5.1f\n", level, iter+1, mse_f, mse_q, mse_m, mse_c);

        //convert_f32matrix_ui8matrix_sat(F1r[level], 0, h-1, 0, w-1, I1r); // saturation necessaire
        //generate_path_filename_k_ndigit_l_extension(dst_path, "F1r_", level, nb_digit, nb_iter+1, "pgm", complete_filename);
        //SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

        //convert_ui16matrix_ui8matrix_sat(E1r[level], 0, h-1, 0, w-1, I1r); // saturation necessaire
        //generate_path_filename_k_ndigit_l_extension(dst_path, "E1r_", level, nb_digit, iter+1, "pgm", complete_filename);
        //SavePGM_ui8matrix(I1r, 0, h-1, 0, w-1, complete_filename);

        if(level) {

            resample_F32(U_f[level], h, w, U_f[level-1]);
            resample_F32(V_f[level], h, w, V_f[level-1]);

            resample_I32(U_q[level], h, w, U_q[level-1]);
            resample_I32(V_q[level], h, w, V_q[level-1]);

            resample_I32(U_m[level], h, w, U_q[level-1]);
            resample_I32(V_m[level], h, w, V_q[level-1]);

            // -- debug F16: OK
            /*if(level == nb_level - 1) {
             generate_path_filename_k_ndigit_extension(dst_path, "Uup_", level, nb_digit, "txt", complete_filename); write_f32matrix(U[level-1], 0, 2*h-1, 0, 2*w-1, "%8.3f", complete_filename);
             generate_path_filename_k_ndigit_extension(dst_path, "Vup_", level, nb_digit, "txt", complete_filename); write_f32matrix(V[level-1], 0, 2*h-1, 0, 2*w-1, "%8.3f", complete_filename);
             }*/
        }

    } // level
end:
    // free
    //puts("free");
    free_ui8matrix(I0,  i0, i1, j0, j1);
    free_ui8matrix(I1,  i0, i1, j0, j1);
    free_ui8matrix(I1r, i0, i1, j0, j1);

    free_f32pyramid(F0,  nb_level, height, width, border);
    free_f32pyramid(F1,  nb_level, height, width, border);
    free_f32pyramid(F1r, nb_level, height, width, border);

    free_f32pyramid(Fx, nb_level, height, width, border);
    free_f32pyramid(Fy, nb_level, height, width, border);
    free_f32pyramid(Ft, nb_level, height, width, border);

    free_f32pyramid(dU_f, nb_level, height, width, border);
    free_f32pyramid(dV_f, nb_level, height, width, border);

    free_f32pyramid(Um_f, nb_level, height, width, border);
    free_f32pyramid(Vm_f, nb_level, height, width, border);

    free_f32pyramid(U_f, nb_level, height, width, border);
    free_f32pyramid(V_f, nb_level, height, width, border);

    free_f32pyramid(Ut_f, nb_level, height, width, 0);
    free_f32pyramid(Vt_f, nb_level, height, width, 0);

    free_ui16pyramid(E0,  nb_level, height, width, border);
    free_ui16pyramid(E1,  nb_level, height, width, border);

    // Qstd
    free_ui16pyramid(E1r_q, nb_level, height, width, border);
    free_ui16pyramid(E1t_q, nb_level, height, width, border);

    free_si16pyramid(Ex_q, nb_level, height, width, border);
    free_si16pyramid(Ey_q, nb_level, height, width, border);
    free_si16pyramid(Et_q, nb_level, height, width, border);

    free_si32pyramid(dU_q, nb_level, height, width, border);
    free_si32pyramid(dV_q, nb_level, height, width, border);

    free_si32pyramid(Um_q, nb_level, height, width, border);
    free_si32pyramid(Vm_q, nb_level, height, width, border);

    free_si32pyramid(U_q, nb_level, height, width, border);
    free_si32pyramid(V_q, nb_level, height, width, border);

    free_si32pyramid(Ut_q, nb_level, height, width, 0);
    free_si32pyramid(Vt_q, nb_level, height, width, 0);

    // Qmax
    free_ui16pyramid(E1r_m, nb_level, height, width, border);
    free_ui16pyramid(E1t_m, nb_level, height, width, border);

    free_si16pyramid(Ex_m, nb_level, height, width, border);
    free_si16pyramid(Ey_m, nb_level, height, width, border);
    free_si16pyramid(Et_m, nb_level, height, width, border);

    free_si32pyramid(dU_m, nb_level, height, width, border);
    free_si32pyramid(dV_m, nb_level, height, width, border);

    free_si32pyramid(Um_m, nb_level, height, width, border);
    free_si32pyramid(Vm_m, nb_level, height, width, border);

    free_si32pyramid(U_m, nb_level, height, width, border);
    free_si32pyramid(V_m, nb_level, height, width, border);

    free_si32pyramid(Ut_m, nb_level, height, width, 0);
    free_si32pyramid(Vt_m, nb_level, height, width, 0);

    free_f32matrix(W, 0, height-1, 0, width-1);

    return 0;
#endif // #ifdef ENABLE_FIXED_POINT
}
// -----------------------------------------
int HS_test_operator(int argc, char *argv[])
// -----------------------------------------
{
    puts("hello !");

    //HS_test_mono1_F32(argc, argv);
    //HS_test_mono2_F32(argc, argv);
    //HS_test_hierarchique1_F32(argc, argv);
    HS_test_hierarchique1_auto_F32(argc, argv);

    //HS_test_mono1_FQ(argc, argv); // comparaison F32, Qstd, Qmax (version 2018)
    //HS_test_hierarchique1_FQ(argc, argv);

    puts("bye...");
    return 0;
}
