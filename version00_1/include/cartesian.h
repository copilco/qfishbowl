#include <iostream>
#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include <fftw3.h>
#include "constant.h"
#include "omp.h"
#include <string>
using namespace std;

//void placeWF(wavefunction wlarge, wavefunction wsmall);
void set_wavefunction_zeros(wavefunction &w);
void placeWF(wavefunction &wlarge, wavefunction &wsmall);
void place_dense_WF(wavefunction &wundense, wavefunction &wdense,  int skiper1, int skiper2, int skiper3);
void gaussian(wavefunction w);
void gaussianR0(wavefunction w,double x1,double x2,double x3);
void gaussianR1(wavefunction &w, double kx1,double kx2,double kx3, double x01,double x02,double x03 );
void gaussian1_wavepacket1D(wavefunction w, double m0, double rho01, double x01, double v01,double b0, double t);
void gaussian_wavepacket1D(wavefunction w, double m0, double rho01, double x01, double v01, double t);
void dipole_XQ1(wavefunction &w, wavefunction &dipolew);
void IMask2D( wavefunction &wfull, wavefunction &wmask,  double _r0, double _r1);


void prop_kinetic(wavefunction &w, hamiltonian &h, complex dt);
void prop_kinetic_complex(wavefunction &w, hamiltonian &h, double dt);
void prop_kinetic_change(wavefunction &w, hamiltonian &h, double dt);
void prop_kinetic_laser_AP(wavefunction &w, hamiltonian &h, complex dt, double vect_pot1, double vect_pot2, double vect_pot3);
void prop_kinetic_laser_AP_change(wavefunction &w, hamiltonian &h, double dt, double vect_pot1, double vect_pot2, double vect_pot3);
void prop_potential(wavefunction &w, hamiltonian h, complex dt);
void prop_potential_complex(wavefunction &w, hamiltonian h, double dt);
void prop_potential_change(wavefunction &w, hamiltonian h, double dt );
void prop_potential_length_gauge(wavefunction &w, hamiltonian h, double dt, double field_e1, double field_e2, double field_e3 ) ;//For propagation in real tim


double kinetic_finite_diff(wavefunction w, hamiltonian h);
double potential_energy(wavefunction w, hamiltonian h);
double potential_energy_length_gauge1D(wavefunction w, hamiltonian h,double efield1);

void momentum_distribution(wavefunction &w , wavefunction &momentum_w , double x0, double y0, double z0);
void momentum_distribution1D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1);
void momentum_distribution2D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1 );
void mask_function1D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1, double _sigma );

void Mask1D( wavefunction &w1 ,wavefunction &momentum_w, double _r0, double _r1);
void IMask1D( wavefunction &wmask, wavefunction &wfull, double _r0, double _r1);
void Mask_Shift1D( wavefunction &w1, wavefunction &momentum_w, double _r0, double _r1, double R0 );
void Mask2D( wavefunction &w1 ,wavefunction &momentum_w, double _r0, double _r1 );
void AntiMask2D( wavefunction &w1 , double _r0, double _r1, double _sigma0 );
double wlecture(FILE *file, wavefunction &w);

double q_expected_kinetic(wavefunction &w, hamiltonian &h);
double q_expected_kinetic_copy(wavefunction &w, hamiltonian &h);
double ionization_scatteringW1D(complex *proj,  int Nmom, double dk);
double position_population(wavefunction &w, double *x, double *y, double *z);
void absorber(wavefunction w, double _frac_x1_left,double _frac_x2_left,double _frac_x3_left,double _frac_x1_right,double _frac_x2_right, double _frac_x3_right, double _exponent);



void linear_combination(wavefunction &w1, wavefunction &w2, wavefunction &w3, vector<double> &b );
void complex_linear_combination(wavefunction &w1, wavefunction &w2, wavefunction &w3, vector<complex> &b );


void project_out(wavefunction &w1, wavefunction &w2);
double bound_population1D(wavefunction &w, wavefunction *warray, int Nstates, int index_initial_state);

complex scattering_wave_projector1D(wavefunction &w, hamiltonian &h, double kx );
void set_scattering_wave_projector1D(wavefunction &w, hamiltonian &h, complex *proj, int Nmom, double dk );
complex projection(wavefunction &w2, wavefunction &w1);
complex projection_on_range(wavefunction &w2, wavefunction &w1,double *x,double *y,double *z);
complex def_integrate( wavefunction &w1, double a1, double b1, double a2, double b2, double a3, double b3);



void q_order(grid &g);
void w_order1D(wavefunction &w, int sw0);
void w_order2D(wavefunction &w, int sw0);

void q_disorder(grid &g);
void w_disorder1D(wavefunction &w, int sw0);
void w_disorder2D(wavefunction &w, int sw0);


//void wavecopy(wavefunction &w1, wavefunction &w2);
void wave_phase2(wavefunction &w);
void snapshot(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3);
void complex_snapshot_WF(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3);


void Qsnapshot(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3);
void Qsnapshot_WF(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3);
void complex_Qsnapshot_WF(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3);


void complex_snapshot_scatteringW1D(FILE *file, wavefunction &w, hamiltonian &h, complex *proj, int Nmom, double dk);
void esay_complex_snapshot_scatteringW1D(FILE *file, wavefunction &w, hamiltonian &h, int Nmom, double dk);


void snapshot_scatteringW1D(FILE *file, wavefunction &w, hamiltonian &h, complex *proj, int Nmom, double dk);
void easy_snapshot_scatteringW1D(FILE *file, wavefunction &w, hamiltonian &h,  int Nmom, double dk);


void snapshot_bound_to_bound_transition(FILE *file, wavefunction &w, wavefunction *warray, int Nstates);
void XQsnapshot(FILE *file, grid &g1, int skiper1, int skiper2, int skiper3);


void biwrite(FILE *file, wavefunction &w);
void biread(FILE *file, wavefunction &w); 
void Rot90(vector<double>& v );


void solver_static_schrodinger_eq(wavefunction &w, hamiltonian &h, double k );
void Continuum_WF( wavefunction &w, hamiltonian &h, double kx );
inline double diff_th( double yn1, double yn_1, double h);
inline double diff(double yn, double yn_1, double h);
inline void cproduct( wavefunction &w, complex a);
inline complex norm_fac( complex psi, complex Dpsi, double k, double xn);



void Coulomb_WF(wavefunction &w, hamiltonian &h, double kx);
double trapz_1D(vector <double> &x, vector<double> &y, double a, double b);


double qmajor(vector<double>& v);
double qminus(vector<double>& v);
double wmax_x(wavefunction &w);
double wmax_q(wavefunction &w);
double hmax_x1D(vector<double>& x, vector<double>& v);
double FWHM(vector<double>& x, vector<double>& v);

double length_au_SI (double _leng0);
double frequency_length_SI (double _w0);
double frequency_au_SI (double _w0);  
double time_au_SI (double _time0);
double intensity_au_SI (double _I0);
double Energy_au_SI(double _E0);

double length_SI_au (double _leng0);
double length_frequency_au (double _lambda);
double frequency_SI_au (double _w0); 
double time_SI_au (double _time0);
double intensity_SI_au (double _I0);

