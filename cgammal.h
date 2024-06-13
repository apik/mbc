#ifndef __CGAMMAL_H__
#define __CGAMMAL_H__
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>
#include "flint/acb.h"

long double arb2ld(const arf_t x, arf_rnd_t rnd)
{
  arf_t t;
  mp_srcptr tp;
  mp_size_t tn;
  long double v;

  arf_init(t);
  /* printf("DOUBLE MAX = %d", __DBL_MANT_DIG__); */
  /* printf("LONG   MAX = %d", __LDBL_MANT_DIG__); */
  /* 53 double, 64 long double */
  arf_set_round(t, x, __LDBL_MANT_DIG__ , rnd);
  ARF_GET_MPN_READONLY(tp, tn, t);
  printf("tn = %d\n", tn);
  if (tn == 1)
    v = (long double)(tp[0]);
  else
    v = (long double)(tp[1]) + (long double)(tp[0]) * ldexpl(1,-64);

  v = ldexpl(v, ARF_EXP(t) - FLINT_BITS);

  if (ARF_SGNBIT(t))
    v = -v;

  arf_clear(t);

  return v;
}


long double complex wgammal(long double complex Z)
{
  /* printf("gamma, z = %LF, %LF\n", creall(Z), cimagl(Z)); */
  /* printf("DOUBLE MAX = %d\n", __DBL_MANT_DIG__); */
  /* printf("LONG   MAX = %d\n", __LDBL_MANT_DIG__); */

  mpfr_t reZmp,imZmp;

  /* 53 double, 64 long double */
  mpfr_init2 (reZmp, 2*__LDBL_MANT_DIG__);
  mpfr_init2 (imZmp, 2*__LDBL_MANT_DIG__);

  /* mpfr_init(reZmp); */
  /* mpfr_init(imZmp); */


  mpfr_set_ld (reZmp, creall(Z), MPFR_RNDZ);
  mpfr_set_ld (imZmp, cimagl(Z), MPFR_RNDZ);

  arf_t arbRe, arbIm;
  arf_init(arbRe);
  arf_init(arbIm);

  arf_set_mpfr(arbRe, reZmp);
  arf_set_mpfr(arbIm, imZmp);

  /* printf("Z="); */
  /* arf_print(arbRe); */
  /* arf_print(arbIm); */
  /* printf("\n"); */


  acb_t arbZ,res;
  acb_init(res);
  acb_init(arbZ);
  /* acb_set_arb_arb(arbZ, arbRe, arbIm); */

  arf_set(arb_midref(acb_realref(arbZ)), arbRe);
  mag_zero(arb_radref(acb_realref(arbZ)));

  arf_set(arb_midref(acb_imagref(arbZ)), arbIm);
  mag_zero(arb_radref(acb_imagref(arbZ)));
  /* printf("Z="); */
  /* acb_print(arbZ); */
  /* printf("\n"); */
  acb_gamma(res, arbZ, 2*__LDBL_MANT_DIG__);
  /* printf("res="); */
  /* acb_print(res); */
  /* printf("\n"); */

  mpfr_t reResmp,imResmp;
  /* 53 double, 64 long double */
  mpfr_init2 (reResmp, 2*__LDBL_MANT_DIG__);
  mpfr_init2 (imResmp, 2*__LDBL_MANT_DIG__);

  arf_get_mpfr(reResmp, arb_midref(acb_realref(res)), MPFR_RNDZ);
  arf_get_mpfr(imResmp, arb_midref(acb_imagref(res)), MPFR_RNDZ);

  long double complex ldRes = mpfr_get_ld(reResmp, MPFR_RNDZ) + I*mpfr_get_ld(imResmp, MPFR_RNDZ);

  acb_clear(res);
  acb_clear(arbZ);

  arf_clear(arbRe);
  arf_clear(arbIm);

  mpfr_clear(reZmp);
  mpfr_clear(imZmp);

  mpfr_clear(reResmp);
  mpfr_clear(imResmp);

  return ldRes;
}


long double complex wpsipgl(long double complex Z, long K)
{

  mpfr_t reZmp,imZmp;
  /* 53 double, 64 long double */
  mpfr_init2 (reZmp, __LDBL_MANT_DIG__);
  mpfr_init2 (imZmp, __LDBL_MANT_DIG__);

  mpfr_set_ld (reZmp, creal(Z), MPFR_RNDZ);
  mpfr_set_ld (imZmp, cimag(Z), MPFR_RNDZ);

  arf_t arbRe, arbIm;
  arf_init(arbRe);
  arf_init(arbIm);

  arf_set_mpfr(arbRe, reZmp);
  arf_set_mpfr(arbIm, imZmp);

  acb_t arbZ,res;
  acb_init(res);
  acb_init(arbZ);

  arf_set(arb_midref(acb_realref(arbZ)), arbRe);
  arf_set(arb_midref(acb_imagref(arbZ)), arbIm);

  acb_t i;
  acb_init(i);

  acb_set_si(i, K);
  acb_polygamma(res, i, arbZ, __LDBL_MANT_DIG__);

  arf_get_mpfr(reZmp, arb_midref(acb_realref(res)), MPFR_RNDZ);
  arf_get_mpfr(imZmp, arb_midref(acb_imagref(res)), MPFR_RNDZ);

  long double complex ldRes = mpfr_get_ld(reZmp, MPFR_RNDZ) + I*mpfr_get_ld(imZmp, MPFR_RNDZ);

  acb_clear(i);
  acb_clear(res);
  acb_clear(arbZ);

  arf_clear(arbRe);
  arf_clear(arbIm);
  mpfr_clear(reZmp);
  mpfr_clear(imZmp);

  return ldRes;
}

#endif  /* __CGAMMA_H__ */
