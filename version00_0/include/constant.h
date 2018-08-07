//Constants
#ifndef CONSTANT_H
#define CONSTANT_H
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <complex>
#define complex complex<double>


const std::complex I=std::complex(0.,1.);
const double lightC_au = 137.036;
const double one_by_lightC_au = 1./lightC_au;

const double pi = 3.141592653589793238463383;
const double dospi = 6.2831853071795862;
const double charge_e1ectron_au = -1.;
const double mass_proton_au=1836.;

enum gauge_t { lengthgauge, velocitygauge, othergauge };

#endif  /* CONSTANT_H */
