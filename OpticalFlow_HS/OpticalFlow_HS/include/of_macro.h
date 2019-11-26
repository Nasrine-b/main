/* ================= */
/* === OF_macro.h === */
/* ================= */

#ifndef __OF_MACRO_H__
#define __OF_MACRO_H__

#define HS_SIZE_MIN 32
#define HS_SIZE_MAX (2*1024)
#define HS_SIZE_STEP 4

#define RUN 2
#define ITER 4


#define RUN 2
#define ITER 4

#define ENABLE_CHRONO_CYCLE
//#define ENABLE_CHRONO_SECONDE

#ifdef ENABLE_CHRONO_CYCLE
#define DUP10(X) X; X; X; X; X; X; X; X; X; X
#define DUP20(X) DUP10(X); DUP10(X)
#define BENCH(X, cpp) tmin=1e30; for(r=0; r<run; r++) { t0=_rdtsc(); DUP20(X); t1=_rdtsc(); dt=(t1-t0)/20.0; if(dt<tmin) tmin = dt;} cpp = tmin/(size*size);
#endif

#ifdef ENABLE_CHRONO_SECONDE
#define BENCH(X, cpp) tmin=1e20; for(r=0; r<run; r++) { t0=dtime();  for(k=0; k<iter; k++) { X; } t1=dtime(); dt=t1-t0; dt /= (double)iter; if(dt<tmin) tmin = dt;} cpp = tmin*FREQ/(size*size);
#endif


#define SWAP_F32(x, y)  do{ float32** t = y; y = x; x = t; }while(0)

#define SWAP_UI8(x, y)  do{  uint8**  t = y; y = x; x = t; }while(0)
#define SWAP_UI16(x, y) do{ uint16**  t = y; y = x; x = t; }while(0)
#define SWAP_SI16(x, y) do{ sint16**  t = y; y = x; x = t; }while(0)
#define SWAP_SI32(x, y) do{ sint32**  t = y; y = x; x = t; }while(0)
#define SWAP_SI64(x, y) do{ sint64**  t = y; y = x; x = t; }while(0)

#define CALC_Q(q) (1 << (q))
//#define CALC_Q(q) pow(2, q)
//#define CALC_Q(q) pow(10, q)

#endif // __OF_MACRO_H__