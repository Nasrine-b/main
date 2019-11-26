/* ================== */
/* === OF_score.c === */
/* ================== */

/*
 * Calcul et affichage des scores et des mesures d'erreur
 */

#ifndef __OF_SCORE_H__
#define __OF_SCORE_H__

#ifdef __cplusplus
extern "C" {
#endif

void OF_calc_scoreRec(OF *of);
void OF_calc_sequenceScoreRec(OF *of);

void OF_display_scoreRec(OF *of);
void OF_display_scoreRecAll(OF *of);
void OF_display_sequenceScoreRec(OF *of);
void OF_display_sequenceScoreRecAll(OF *of);

#ifdef __cplusplus
}
#endif

#endif // __OF_SCORE_H__