#include <iostream>

void placeWF(wavefunction wlarge, wavefunction wsmall);
void gaussian(wavefunction w);
void gaussian(wavefunction w,double x1,double x2,double x3);
void gaussianR1(wavefunction w, double kx1,double kx2,double kx3, double x01,double x02,double x03 );

void prop_kinetic(wavefunction &w, hamiltonian &h, complex dt);
void prop_kinetic_laser_AP(wavefunction &w, hamiltonian &h, complex dt, double vect_pot1, double vect_pot2, double vect_pot3);
void prop_potential(wavefunction &w, hamiltonian h, complex dt);
void prop_potential_laser_ER(wavefunction &w, hamiltonian h, complex dt, double e_field1, double e_field2, double e_field3);

double kinetic_finite_diff(wavefunction w, hamiltonian h);
double potential_energy(wavefunction w, hamiltonian h);

void momentum_distribution(wavefunction &w , wavefunction &momentum_w , double x0, double y0, double z0);
void momentum_distribution1D(wavefunction &w , wavefunction &momentum_w ,  double _r0, double _r1);
void mask_function1D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1, double _sigma );

void Mask2D( wavefunction &w1 , double _r0, double _r1, double _sigma0 );
void AntiMask2D( wavefunction &w1 , double _r0, double _r1, double _sigma0 );

double q_expected_kinetic(wavefunction &w, hamiltonian &h);
void absorber(wavefunction w, double _frac_x1_left,double _frac_x2_left,double _frac_x3_left,double _frac_x1_right,double _frac_x2_right, double _frac_x3_right, double _exponent);

complex projection(wavefunction &w1, wavefunction &w2);
void project_out(wavefunction &w1, wavefunction &w2);

void snapshot(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3);
void Qsnapshot(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3);
