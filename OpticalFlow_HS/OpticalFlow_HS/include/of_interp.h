/* =================== */
/* === OF_interp.h === */
/* =================== */

#ifndef __OF_INTERP_H__
#define __OF_INTERP_H__

#ifdef __cplusplus
extern "C" {
#endif
    
void FOLKI_interpLinear   (OF *of);
void FOLKI_interpCubic    (OF *of);

void HS_interpLinear   (OF *of);
void HS_interpCubic    (OF *of);

void FOLKI_interpLinearLUT(OF *of);
void FOLKI_interpCubicLUT (OF *of);
    
void FOLKI_gradientSobel(OF *of);

void FOLKI_gradientLinear   (OF *of);
void FOLKI_gradientCubic    (OF *of);
void FOLKI_gradientLinearLUT(OF *of);
void FOLKI_gradientCubicLUT (OF *of);
    
#ifdef __cplusplus
}
#endif

#endif // __OF_INTERP_H__