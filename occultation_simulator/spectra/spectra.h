#ifndef _SPECTRA_H_
#define _SPECTRA_H_

#define ABSMAG -0.03579         // 0.264
#define RADIUS 1.724            // 0.036

#include "v_filter.h"

#include "a0i.h"
#include "a0iii.h"
#include "a0iv.h"
#include "a0v.h"
#include "a2i.h"
#include "a2v.h"
#include "a3iii.h"
#include "a3v.h"
#include "a47iv.h"
#include "a5iii.h"
#include "a5v.h"
#include "a7iii.h"
#include "a7v.h"
#include "b0i.h"
#include "b0v.h"
#include "b12iii.h"
#include "b1i.h"
#include "b1v.h"
#include "b2ii.h"
#include "b2iv.h"
#include "b3i.h"
#include "b3iii.h"
#include "b3v.h"
#include "b57v.h"
#include "b5i.h"
#include "b5ii.h"
#include "b5iii.h"
#include "b6iv.h"
#include "b8i.h"
#include "b8v.h"
#include "b9iii.h"
#include "b9v.h"
#include "f02iv.h"
#include "f0i.h"
#include "f0ii.h"
#include "f0iii.h"
#include "f0v.h"
#include "f2ii.h"
#include "f2iii.h"
#include "f2v.h"
#include "f5i.h"
#include "f5iii.h"
#include "f5iv.h"
#include "f5v.h"
#include "f6v.h"
#include "f8i.h"
#include "f8iv.h"
#include "f8v.h"
#include "g0i.h"
#include "g0iii.h"
#include "g0iv.h"
#include "g0v.h"
#include "g2i.h"
#include "g2iv.h"
#include "g2v.h"
#include "g5i.h"
#include "g5ii.h"
#include "g5iii.h"
#include "g5iv.h"
#include "g5v.h"
#include "g8i.h"
#include "g8iii.h"
#include "g8iv.h"
#include "g8v.h"
#include "k01ii.h"
#include "k0iii.h"
#include "k0iv.h"
#include "k0v.h"
#include "k1iii.h"
#include "k1iv.h"
#include "k2i.h"
#include "k2iii.h"
#include "k2v.h"
#include "k34ii.h"
#include "k3i.h"
#include "k3iii.h"
#include "k3iv.h"
#include "k3v.h"
#include "k4i.h"
#include "k4iii.h"
#include "k4v.h"
#include "k5iii.h"
#include "k5v.h"
#include "k7v.h"
#include "m0iii.h"
#include "m0v.h"
#include "m10iii.h"
#include "m1iii.h"
#include "m1v.h"
#include "m2.5v.h"
#include "m2i.h"
#include "m2iii.h"
#include "m2v.h"
#include "m3ii.h"
#include "m3iii.h"
#include "m3v.h"
#include "m4iii.h"
#include "m4v.h"
#include "m5iii.h"
#include "m5v.h"
#include "m6iii.h"
#include "m6v.h"
#include "m7iii.h"
#include "m8iii.h"
#include "m9iii.h"
#include "o5v.h"
#include "o8iii.h"
#include "o9v.h"
#include "fltopt.h"
#include "ScoX1.h"

static double *_lambda[] = {
   _o5v_lambda,
   _o8iii_lambda,
   _o9v_lambda,
   _b0i_lambda,
   _b0v_lambda,
   _b12iii_lambda,
   _b1i_lambda,
   _b1v_lambda,
   _b2ii_lambda,
   _b2iv_lambda,
   _b3i_lambda,
   _b3iii_lambda,
   _b3v_lambda,
   _b57v_lambda,
   _b5i_lambda,
   _b5ii_lambda,
   _b5iii_lambda,
   _b6iv_lambda,
   _b8i_lambda,
   _b8v_lambda,
   _b9iii_lambda,
   _b9v_lambda,
   _a0i_lambda,
   _a0iii_lambda,
   _a0iv_lambda,
   _a0v_lambda,
   _a2i_lambda,
   _a2v_lambda,
   _a3iii_lambda,
   _a3v_lambda,
   _a47iv_lambda,
   _a5iii_lambda,
   _a5v_lambda,
   _a7iii_lambda,
   _a7v_lambda,
   _f02iv_lambda,
   _f0i_lambda,
   _f0ii_lambda,
   _f0iii_lambda,
   _f0v_lambda,
   _f2ii_lambda,
   _f2iii_lambda,
   _f2v_lambda,
   _f5i_lambda,
   _f5iii_lambda,
   _f5iv_lambda,
   _f5v_lambda,
   _f6v_lambda,
   _f8i_lambda,
   _f8iv_lambda,
   _f8v_lambda,
   _g0i_lambda,
   _g0iii_lambda,
   _g0iv_lambda,
   _g0v_lambda,
   _g2i_lambda,
   _g2iv_lambda,
   _g2v_lambda,
   _g5i_lambda,
   _g5ii_lambda,
   _g5iii_lambda,
   _g5iv_lambda,
   _g5v_lambda,
   _g8i_lambda,
   _g8iii_lambda,
   _g8iv_lambda,
   _g8v_lambda,
   _k01ii_lambda,
   _k0iii_lambda,
   _k0iv_lambda,
   _k0v_lambda,
   _k1iii_lambda,
   _k1iv_lambda,
   _k2i_lambda,
   _k2iii_lambda,
   _k2v_lambda,
   _k34ii_lambda,
   _k3i_lambda,
   _k3iii_lambda,
   _k3iv_lambda,
   _k3v_lambda,
   _k4i_lambda,
   _k4iii_lambda,
   _k4v_lambda,
   _k5iii_lambda,
   _k5v_lambda,
   _k7v_lambda,
   _m0iii_lambda,
   _m0v_lambda,
   _m10iii_lambda,
   _m1iii_lambda,
   _m1v_lambda,
   _m2_5v_lambda,
   _m2i_lambda,
   _m2iii_lambda,
   _m2v_lambda,
   _m3ii_lambda,
   _m3iii_lambda,
   _m3v_lambda,
   _m4iii_lambda,
   _m4v_lambda,
   _m5iii_lambda,
   _m5v_lambda,
   _m6iii_lambda,
   _m6v_lambda,
   _m7iii_lambda,
   _m8iii_lambda,
   _m9iii_lambda,
   _fltopt_lambda,
   _scox1_lambda
};

static double *_flux[] = {
   _o5v_flux,
   _o8iii_flux,
   _o9v_flux,
   _b0i_flux,
   _b0v_flux,
   _b12iii_flux,
   _b1i_flux,
   _b1v_flux,
   _b2ii_flux,
   _b2iv_flux,
   _b3i_flux,
   _b3iii_flux,
   _b3v_flux,
   _b57v_flux,
   _b5i_flux,
   _b5ii_flux,
   _b5iii_flux,
   _b6iv_flux,
   _b8i_flux,
   _b8v_flux,
   _b9iii_flux,
   _b9v_flux,
   _a0i_flux,
   _a0iii_flux,
   _a0iv_flux,
   _a0v_flux,
   _a2i_flux,
   _a2v_flux,
   _a3iii_flux,
   _a3v_flux,
   _a47iv_flux,
   _a5iii_flux,
   _a5v_flux,
   _a7iii_flux,
   _a7v_flux,
   _f02iv_flux,
   _f0i_flux,
   _f0ii_flux,
   _f0iii_flux,
   _f0v_flux,
   _f2ii_flux,
   _f2iii_flux,
   _f2v_flux,
   _f5i_flux,
   _f5iii_flux,
   _f5iv_flux,
   _f5v_flux,
   _f6v_flux,
   _f8i_flux,
   _f8iv_flux,
   _f8v_flux,
   _g0i_flux,
   _g0iii_flux,
   _g0iv_flux,
   _g0v_flux,
   _g2i_flux,
   _g2iv_flux,
   _g2v_flux,
   _g5i_flux,
   _g5ii_flux,
   _g5iii_flux,
   _g5iv_flux,
   _g5v_flux,
   _g8i_flux,
   _g8iii_flux,
   _g8iv_flux,
   _g8v_flux,
   _k01ii_flux,
   _k0iii_flux,
   _k0iv_flux,
   _k0v_flux,
   _k1iii_flux,
   _k1iv_flux,
   _k2i_flux,
   _k2iii_flux,
   _k2v_flux,
   _k34ii_flux,
   _k3i_flux,
   _k3iii_flux,
   _k3iv_flux,
   _k3v_flux,
   _k4i_flux,
   _k4iii_flux,
   _k4v_flux,
   _k5iii_flux,
   _k5v_flux,
   _k7v_flux,
   _m0iii_flux,
   _m0v_flux,
   _m10iii_flux,
   _m1iii_flux,
   _m1v_flux,
   _m2_5v_flux,
   _m2i_flux,
   _m2iii_flux,
   _m2v_flux,
   _m3ii_flux,
   _m3iii_flux,
   _m3v_flux,
   _m4iii_flux,
   _m4v_flux,
   _m5iii_flux,
   _m5v_flux,
   _m6iii_flux,
   _m6v_flux,
   _m7iii_flux,
   _m8iii_flux,
   _m9iii_flux,
   _fltopt_flux,
   _scox1_flux
};

/*enum, id, filterstring, abs mag, default app mag, effective temp in K, star radius in R_sun, ang size (in rad?), flux entries, &wavelength, &flux*/

static struct spectype _spectypeinfo[] = {
   {O5V, "o5v", "", -5.7, 12, 42000, 12, -1, 0.0, 1895, NULL, NULL},
   {O8III, "o8iii", "", -5.8, 12, 34700, 30, -1, 0.0, 1895, NULL, NULL},
   {O9V, "o9v", "", -4.5, 12, 34000, -1, -1, 0.0, 1895, NULL, NULL},
   {B0I, "b0i", "", -6.4, 12, 26000, 30, -1, 0.0, 1895, NULL, NULL},
   {B0V, "b0v", "", -4, 12, 30000, 7.4, -1, 0.0, 1895, NULL, NULL},
   {B12III, "b12iii", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {B1I, "b1i", "", -6.4, 12, 20800, -1, -1, 0.0, 1895, NULL, NULL},
   {B1V, "b1v", "", -3.2, 12, 25400, -1, -1, 0.0, 1895, NULL, NULL},
   {B2II, "b2ii", "", -4.8, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {B2IV, "b2iv", "", -3.1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {B3I, "b3i", "", -6.3, 12, 16200, -1, -1, 0.0, 1895, NULL, NULL},
   {B3III, "b3iii", "", -3, 12, 17100, -1, -1, 0.0, 1895, NULL, NULL},
   {B3V, "b3v", "", -1.6, 12, 18700, 4.8, -1, 0.0, 1895, NULL, NULL},
   {B57V, "b57v", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {B5I, "b5i", "", -6.2, 12, 13600, 50, -1, 0.0, 1895, NULL, NULL},
   {B5II, "b5ii", "", -4, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {B5III, "b5iii", "", -2.2, 12, 15000, 8, -1, 0.0, 1895, NULL, NULL},
   {B6IV, "b6iv", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {B8I, "b8i", "", -6.2, 12, 11200, -1, -1, 0.0, 1895, NULL, NULL},
   {B8V, "b8v", "", -0.25, 12, 11400, 3, -1, 0.0, 1895, NULL, NULL},
   {B9III, "b9iii", "", -0.6, 12, 11000, -1, -1, 0.0, 1895, NULL, NULL},
   {B9V, "b9v", "", 0.2, 12, 10500, -1, -1, 0.0, 1895, NULL, NULL},
   {A0I, "a0i", "", -6.3, 12, 9980, 60, -1, 0.0, 1895, NULL, NULL},
   {A0III, "a0iii", "", 0, 12, 10100, 5, -1, 0.0, 1895, NULL, NULL},
   {A0IV, "a0iv", "", 0.3, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {A0V, "a0v", "", 0.65, 12, 9790, 2.4, -1, 0.0, 1895, NULL, NULL},
   {A2I, "a2i", "", -6.5, 12, 9080, -1, -1, 0.0, 1895, NULL, NULL},
   {A2V, "a2v", "", 1.3, 12, 9000, -1, -1, 0.0, 1895, NULL, NULL},
   {A3III, "a3iii", "", 0.5, 12, 8600, -1, -1, 0.0, 1895, NULL, NULL},
   {A3V, "a3v", "", 1.5, 12, 8720, -1, -1, 0.0, 1895, NULL, NULL},
   {A47IV, "a47iv", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {A5III, "a5iii", "", 0.7, 12, 8100, -1, -1, 0.0, 1895, NULL, NULL},
   {A5V, "a5v", "", 1.95, 12, 8180, 1.7, -1, 0.0, 1895, NULL, NULL},
   {A7III, "a7iii", "", 1.1, 12, 7650, -1, -1, 0.0, 1895, NULL, NULL},
   {A7V, "a7v", "", 2.2, 12, 7850, -1, -1, 0.0, 1895, NULL, NULL},
   {F02IV, "f02iv", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {F0I, "f0i", "", -6.6, 12, 7460, 80, -1, 0.0, 1895, NULL, NULL},
   {F0II, "f0ii", "", -2.5, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {F0III, "f0iii", "", 1.5, 12, 7150, -1, -1, 0.0, 1895, NULL, NULL},
   {F0V, "f0v", "", 2.7, 12, 7300, 1.5, -1, 0.0, 1895, NULL, NULL},
   {F2II, "f2ii", "", -2.4, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {F2III, "f2iii", "", 1.7, 12, 6870, -1, -1, 0.0, 1895, NULL, NULL},
   {F2V, "f2v", "", 3.6, 12, 7000, -1, -1, 0.0, 1895, NULL, NULL},
   {F5I, "f5i", "", -6.6, 12, 6370, 100, -1, 0.0, 1895, NULL, NULL},
   {F5III, "f5iii", "", 1.6, 12, 6470, -1, -1, 0.0, 1895, NULL, NULL},
   {F5IV, "f5iv", "", 2.5, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {F5V, "f5v", "", 3.5, 12, 6650, 1.3, -1, 0.0, 1895, NULL, NULL},
   {F6V, "f6v", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {F8I, "f8i", "", -6.5, 12, 5750, -1, -1, 0.0, 1895, NULL, NULL},
   {F8IV, "f8iv", "", 2.8, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {F8V, "f8v", "", 4, 12, 6250, -1, -1, 0.0, 1895, NULL, NULL},
   {G0I, "g0i", "", -6.4, 12, 5370, 120, -1, 0.0, 1895, NULL, NULL},
   {G0III, "g0iii", "", 1, 12, 5850, 6, -1, 0.0, 1895, NULL, NULL},
   {G0IV, "g0iv", "", 3, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {G0V, "g0v", "", 4.4, 12, 5940, 1.1, -1, 0.0, 1895, NULL, NULL},
   {G2I, "g2i", "", -6.3, 12, 5190, -1, -1, 0.0, 1895, NULL, NULL},
   {G2IV, "g2iv", "", 3, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {G2V, "g2v", "", 4.7, 12, 5790, -1, -1, 0.0, 1895, NULL, NULL},
   {G5I, "g5i", "", -6.2, 12, 4930, 150, -1, 0.0, 1895, NULL, NULL},
   {G5II, "g5ii", "", -2.3, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {G5III, "g5iii", "", 0.9, 12, 5050, 10, -1, 0.0, 1895, NULL, NULL},
   {G5IV, "g5iv", "", 3.1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {G5V, "g5v", "", 5.1, 12, 5560, 0.92, -1, 0.0, 1895, NULL, NULL},
   {G8I, "g8i", "", -6.1, 12, 4700, -1, -1, 0.0, 1895, NULL, NULL},
   {G8III, "g8iii", "", 0.8, 12, 4800, -1, -1, 0.0, 1895, NULL, NULL},
   {G8IV, "g8iv", "", 3.1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {G8V, "g8v", "", 5.5, 12, 5310, -1, -1, 0.0, 1895, NULL, NULL},
   {K01II, "k01ii", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {K0III, "k0iii", "", 0.7, 12, 4660, 15, -1, 0.0, 1895, NULL, NULL},
   {K0IV, "k0iv", "", 3.1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {K0V, "k0v", "", 5.9, 12, 5150, 0.85, -1, 0.0, 1895, NULL, NULL},
   {K1III, "k1iii", "", 0.6, 12, 4600, -1, -1, 0.0, 1895, NULL, NULL},
   {K1IV, "k1iv", "", 3.1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {K2I, "k2i", "", -5.9, 12, 4310, -1, -1, 0.0, 1895, NULL, NULL},
   {K2III, "k2iii", "", 0.5, 12, 4390, -1, -1, 0.0, 1895, NULL, NULL},
   {K2V, "k2v", "", 6.4, 12, 4830, -1, -1, 0.0, 1895, NULL, NULL},
   {K34II, "k34ii", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {K3I, "k3i", "", -5.9, 12, 4080, -1, -1, 0.0, 1895, NULL, NULL},
   {K3III, "k3iii", "", 0.3, 12, 4200, -1, -1, 0.0, 1895, NULL, NULL},
   {K3IV, "k3iv", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {K3V, "k3v", "", 6.65, 12, 4730, -1, -1, 0.0, 1895, NULL, NULL},
   {K4I, "k4i", "", -5.8, 12, 3950, -1, -1, 0.0, 1895, NULL, NULL},
   {K4III, "k4iii", "", 0, 12, 4000, -1, -1, 0.0, 1895, NULL, NULL},
   {K4V, "k4v", "", 7, 12, 4590, -1, -1, 0.0, 1895, NULL, NULL},
   {K5III, "k5iii", "", -0.2, 12, 4050, 25, -1, 0.0, 1895, NULL, NULL},
   {K5V, "k5v", "", 7.35, 12, 4410, 0.72, -1, 0.0, 1895, NULL, NULL},
   {K7V, "k7v", "", 8.1, 12, 4060, -1, -1, 0.0, 1895, NULL, NULL},
   {M0III, "m0iii", "", -0.4, 12, 3690, 40, -1, 0.0, 1895, NULL, NULL},
   {M0V, "m0v", "", 8.8, 12, 3840, 0.6, -1, 0.0, 1895, NULL, NULL},
   {M10III, "m10iii", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {M1III, "m1iii", "", -1, 12, 3720, -1, -1, 0.0, 1895, NULL, NULL},
   {M1V, "m1v", "", 9.3, 12, 3720, -1, -1, 0.0, 1895, NULL, NULL},
   {M2_5V, "m2.5v", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {M2I, "m2i", "", -5.6, 12, 3370, 800, -1, 0.0, 1895, NULL, NULL},
   {M2III, "m2iii", "", -0.6, 12, 3540, -1, -1, 0.0, 1895, NULL, NULL},
   {M2V, "m2v", "", 9.9, 12, 3520, 0.5, -1, 0.0, 1895, NULL, NULL},
   {M3II, "m3ii", "", -2.6, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {M3III, "m3iii", "", -0.6, 12, 3530, -1, -1, 0.0, 1895, NULL, NULL},
   {M3V, "m3v", "", 10.4, 12, 3470, -1, -1, 0.0, 1895, NULL, NULL},
   {M4III, "m4iii", "", -1, 12, 3430, -1, -1, 0.0, 1895, NULL, NULL},
   {M4V, "m4v", "", 11.3, 12, 3370, -1, -1, 0.0, 1895, NULL, NULL},
   {M5III, "m5iii", "", -0.3, 12, 3380, -1, -1, 0.0, 1895, NULL, NULL},
   {M5V, "m5v", "", 12.3, 12, 3170, 0.27, -1, 0.0, 1895, NULL, NULL},
   {M6III, "m6iii", "", -1, 12, 3240, -1, -1, 0.0, 1895, NULL, NULL},
   {M6V, "m6v", "", -1, 12, 3050, -1, -1, 0.0, 1895, NULL, NULL},
   {M7III, "m7iii", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {M8III, "m8iii", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {M9III, "m9iii", "", -1, 12, -1, -1, -1, 0.0, 1895, NULL, NULL},
   {FLTOPT, "fltopt", "", ABSMAG, 12.0, -1, RADIUS, 1.4e-11, 0.0, 1895, NULL,
    NULL},
   {SCOX1, "scox1", "", ABSMAG, 12.0, -1, RADIUS, 1.4e-11, 0.0, 1895, NULL,
    NULL}
};

#endif
