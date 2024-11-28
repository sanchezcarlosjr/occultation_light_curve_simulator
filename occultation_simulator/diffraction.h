//
// Created by cest on 1/31/24.
//

#ifndef OCCULTATION_LIGHT_CURVES_DIFFRACTION_H
#define OCCULTATION_LIGHT_CURVES_DIFFRACTION_H
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

double SNR_TAOS2(double mV);
double calcPlano(double d, double lmda, double ua);
double astronomicalUnitsToMeters(double ua);
gsl_matrix* pupilCO(int M, double D, double d);
gsl_matrix* pupilDoble(int M, double D, double d);
void promedioPD(gsl_matrix_complex *diffractionPattern, double R_star, double plano, int M, double d, gsl_matrix_complex *intensityOut);

#endif //OCCULTATION_LIGHT_CURVES_DIFFRACTION_H
