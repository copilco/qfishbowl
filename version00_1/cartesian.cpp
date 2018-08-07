#include "cartesian.h"


void set_wavefunction_zeros(wavefunction &w){
	
	for(int k=0;k<w.n3*w.n2*w.n1;k++)
	{ 
		w.w[k][0] = w.w[k][1] = 0.;
	}
	
}

//Put wavefunction wsmall on a big wavefunction wlarge
void placeWF(wavefunction &wlarge,  wavefunction &wsmall)
{
	int kplacer=floor((wlarge.n3-wsmall.n3)/2);
	int jplacer=floor((wlarge.n2-wsmall.n2)/2);
	int iplacer=floor((wlarge.n1-wsmall.n1)/2);

	if(wsmall.n1!=1 && iplacer<0){
		cout << "\n\n*** Check wlarge wavefunction size grid *** \n\n";
		exit (1);
	}
	
	set_wavefunction_zeros(wlarge);
	
	for(int k=0;k<wsmall.n3;k++)
		for(int j=0;j<wsmall.n2;j++)
			for(int i=0;i<wsmall.n1;i++)
			{ 	  
				wlarge.w[wlarge.index(kplacer+k,jplacer+j,iplacer+i)][0]
						=	wsmall.w[wsmall.index(k,j,i)][0];
				wlarge.w[wlarge.index(kplacer+k,jplacer+j,iplacer+i)][1]
						=	wsmall.w[wsmall.index(k,j,i)][1];
			}
}



//Gaussian function version 0
void gaussian(wavefunction w)
{
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{ 	  
				w.w[w.index(k,j,i)][0] = exp( -w.x1[i]*w.x1[i]-w.x2[j]*w.x2[j]-w.x3[k]*w.x3[k] );//*sin(2.2*w.x2[j]);
				//*cos(2.2*w.x2[j]);
				w.w[w.index(k,j,i)][1] = 0.;//  exp( -w.x1[i]*w.x1[i]-w.x2[j]*w.x2[j]-w.x3[k]*w.x3[k] )*sin(2.2*w.x2[j]);;
			}
}



//Gaussian function version 1
void gaussianR0(wavefunction w, double x01,double x02,double x03 )
{
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{ 	  
				w.w[w.index(k,j,i)][0] = exp( -(w.x1[i]-x01)*(w.x1[i]-x01) - (w.x2[j]-x02)*(w.x2[j]-x02)
											- (w.x3[k]-x03)*(w.x3[k]-x03) );//*cos(4*w.x1[i] );
				w.w[w.index(k,j,i)][1] = 0.;//exp( -(w.x1[i]-x01)*(w.x1[i]-x01)-(w.x2[j]-x02)*(w.x2[j]-x02)-(w.x3[k]-x03)*(w.x3[k]-x03) )*sin(4*w.x1[j]);
			}
}



//Gaussian function version 2
void gaussianR1(wavefunction &w, double kx1,double kx2,double kx3, double x01,double x02,double x03 )
{
	double rhox1=1.;
	double rhox2=1.;
	double rhox3=1.;
	int sw0=0;
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{ 	  
				w.w[w.index(k,j,i)][0]= exp( - (w.x1[i]-x01)*(w.x1[i]-x01)*rhox1*rhox1 - sw0*(w.x2[j]-x02)*(w.x2[j]-x02)*rhox2*rhox2
											 - sw0*(w.x3[k]-x03)*(w.x3[k]-x03)*rhox3*rhox3 )*cos(kx1*w.x1[i] + sw0*kx2*w.x2[j] + sw0*kx3*w.x3[k] );

				w.w[w.index(k,j,i)][1]= exp( - (w.x1[i]-x01)*(w.x1[i]-x01)*rhox1*rhox1-(w.x2[j]-x02)*(w.x2[j]-x02)*rhox2*rhox2
											 - (w.x3[k]-x03)*(w.x3[k]-x03)*rhox3*rhox3 )*sin(kx1*w.x1[i] + sw0*kx2*w.x2[j] + sw0*kx3*w.x3[k]);
			}
}




//Gaussian function version 3
void gaussian1_wavepacket1D(wavefunction w, double m0, double rho01, double x01, double v01, double b0, double t)
{
	double k01=m0*v01; //abs(m0*v01); //m0*v01;//
	double a=2.*rho01;
	double Ek=1/2.*m0*v01*v01;
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{ 	  
				w.w[i][0]= exp( -(w.x1[i]-v01*t-x01)*(w.x1[i]-v01*t-x01)/a/a )*cos(k01*(w.x1[i]-x01) +b0*(w.x1[i]-x01)*(w.x1[i]-x01)-Ek*t );
				w.w[i][1]= exp( -(w.x1[i]-v01*t-x01)*(w.x1[i]-v01*t-x01)/a/a )*sin(k01*(w.x1[i]-x01) +b0*(w.x1[i]-x01)*(w.x1[i]-x01)-Ek*t );
			}
}




//Gaussian function version 4
void gaussian_wavepacket1D(wavefunction w, double m0, double rho01, double x01, double v01, double t)
{
    double a = 2.*rho01;
    double w0=1.;
    double b= 2./m0;
	//v01=v01*cos( w0*t +pi/4)*cos( w0*t +pi/4);
    double k01=    m0*v01;
    double thetat  = 1./2.*atan2(b*t,a*a);
    double phit    = - thetat - k01*v01*t/2.;
    double gammat  = 1./( pow(a,4.) + pow(b*t,2.) );
	
    double ampt    = pow( (2.*a*a*gammat/pi) , 1./4.) ;    
    for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{   
				double phaset  = b*t*(w.x1[i]-v01*t-x01)*(w.x1[i]-v01*t-x01)*gammat + 
				k01*(w.x1[i]-x01) + phit;
				
				w.w[w.index(k,j,i)][0]= ampt*exp( -a*a*(w.x1[i]-v01*t-x01)*(w.x1[i]-v01*t-x01)*gammat )*cos(phaset);
				w.w[w.index(k,j,i)][1]= ampt*exp( -a*a*(w.x1[i]-v01*t-x01)*(w.x1[i]-v01*t-x01)*gammat )*sin(phaset);
			}
}



//Dipole matrix element by plane wave projection
void dipole_XQ1(wavefunction &w, wavefunction &dipolew)
{
	
    for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				dipolew.w[w.index(k,j,i)][0]= -w.x1[i]*w.w[w.index(k,j,i)][0];
				dipolew.w[w.index(k,j,i)][1]= -w.x1[i]*w.w[w.index(k,j,i)][1];
			}
	dipolew.Direct_Fourier_Transform();
}




//Kinetic Evolver 
void prop_kinetic(wavefunction &w, hamiltonian &h, complex dt)
{
	w.expected_kin=0.;
	w.qnorm=0.;
	double factor=1./w.n1/w.n2/w.n3;
	/**********************/
	fftw_execute(w.p3DF);
	/**********************/  
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{	  
				double kinetic=
				w.q1[i]*w.q1[i]*h.a1+
				w.q2[j]*w.q2[j]*h.a2+
				w.q3[k]*w.q3[k]*h.a3;
				
				complex kinoperator=exp(-I*dt*kinetic );
				
				double aux0r=w.w[w.index(k,j,i)][0];
				double aux0i=w.w[w.index(k,j,i)][1];
				
				w.w[w.index(k,j,i)][0]=( aux0r*real(kinoperator)+aux0i*imag(kinoperator) )*factor;
				w.w[w.index(k,j,i)][1]=( aux0i*real(kinoperator)-aux0r*imag(kinoperator) )*factor;
			}       
	/**********************/
	fftw_execute(w.p3DB);
	/**********************/
}




//Kinetic Evolver Version 2
void prop_kinetic_complex(wavefunction &w, hamiltonian &h, double dt) //For propagation in imaginary time
{
	w.expected_kin=0.;
	w.qnorm=0.;
	double factor=1./w.n1/w.n2/w.n3;
	/**********************/
	fftw_execute(w.p3DF);
	/**********************/  
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				
				double kinetic=
				w.q1[i]*w.q1[i]*h.a1+
				w.q2[j]*w.q2[j]*h.a2+
				w.q3[k]*w.q3[k]*h.a3;
				
				complex kinoperator=complex(exp(-dt*kinetic ),0);
				complex aux=complex(w.w[w.index(k,j,i)][0], w.w[w.index(k,j,i)][1]);
				aux=aux*kinoperator;
				
				w.w[w.index(k,j,i)][0]=real(aux)*factor;
				w.w[w.index(k,j,i)][1]=imag(aux)*factor;
			}       
	/**********************/
	fftw_execute(w.p3DB);
	/**********************/
}




//Kinetic Evolver Version 3
void prop_kinetic_change(wavefunction &w, hamiltonian &h, double dt) //For propagation in real time
{
	w.expected_kin=0.;
	w.qnorm=0.;
	double factor=1./w.n1/w.n2/w.n3;
	/**********************/
	fftw_execute(w.p3DF);
	/**********************/  
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				
				double kinetic=
				w.q1[i]*w.q1[i]*h.a1+
				w.q2[j]*w.q2[j]*h.a2+
				w.q3[k]*w.q3[k]*h.a3;
				
				complex kinoperator=complex(cos(dt*kinetic),-sin(dt*kinetic));
				complex aux=complex( w.w[w.index(k,j,i)][0], w.w[w.index(k,j,i)][1] )*kinoperator;
				
				w.w[w.index(k,j,i)][0] = real(aux)*factor; 
				w.w[w.index(k,j,i)][1] = imag(aux)*factor; 
			}
	
	/**********************/
	fftw_execute(w.p3DB);
	/**********************/
}






//Kinetic Evolver Version 4
void prop_kinetic_laser_AP(wavefunction &w, hamiltonian &h, complex dt, double vect_pot1, double vect_pot2, double vect_pot3) //For propagation in real and imaginary time
{
	
	w.expected_kin=0.;
	w.qnorm=0.;
	
	double factor=1./w.n1/w.n2/w.n3;
	/**********************/
	fftw_execute(w.p3DF);
	/**********************/
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{	  
				double kinetic=
				(w.q1[i]-h.b1*vect_pot1/lightC_au)*(w.q1[i]-h.b1*vect_pot1/lightC_au)*h.a1+
				(w.q2[j]-h.b2*vect_pot2/lightC_au)*(w.q2[j]-h.b2*vect_pot2/lightC_au)*h.a2+
				(w.q3[k]-h.b3*vect_pot3/lightC_au)*(w.q3[k]-h.b3*vect_pot3/lightC_au)*h.a3;
				
				w.expected_kin+=w.dq1*w.dq2*w.dq3*kinetic*
				(w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1])
				*w.kfactor1*w.kfactor1*w.kfactor2*w.kfactor2*w.kfactor3*w.kfactor3;
				
				w.qnorm+=w.dq1*w.dq2*w.dq3*
				(w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1])
				*w.kfactor1*w.kfactor1*w.kfactor2*w.kfactor2*w.kfactor3*w.kfactor3;
				
				complex kinoperator=exp(-I*dt*kinetic);
				
				double aux0r=w.w[w.index(k,j,i)][0];
				double aux0i=w.w[w.index(k,j,i)][1];
				
				w.w[w.index(k,j,i)][0]=( aux0r*real(kinoperator)+aux0i*imag(kinoperator) )*factor;
				w.w[w.index(k,j,i)][1]=( aux0i*real(kinoperator)-aux0r*imag(kinoperator) )*factor;
				
			}
	
	/**********************/
	fftw_execute(w.p3DB);
	/**********************/     
}




//Kinetic Evolver Version 5
void prop_kinetic_laser_AP_change(wavefunction &w, hamiltonian &h, double dt, double vect_pot1, double vect_pot2, double vect_pot3)
{
	
	w.expected_kin=0.;
	w.qnorm=0.;
	
	double factor=1./w.n1/w.n2/w.n3;
	/**********************/
	fftw_execute(w.p3DF);
	/**********************/  
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{	  
				double kinetic=
				(w.q1[i]-h.b1*vect_pot1/lightC_au)*(w.q1[i]-h.b1*vect_pot1/lightC_au)*h.a1+
				(w.q2[j]-h.b2*vect_pot2/lightC_au)*(w.q2[j]-h.b2*vect_pot2/lightC_au)*h.a2+
				(w.q3[k]-h.b3*vect_pot3/lightC_au)*(w.q3[k]-h.b3*vect_pot3/lightC_au)*h.a3;	  
				
				complex kinoperator=complex(cos(dt*kinetic), -sin(dt*kinetic));
				complex aux=complex(w.w[w.index(k,j,i)][0], w.w[w.index(k,j,i)][1])*kinoperator;
				
				w.w[w.index(k,j,i)][0]=real(aux)*factor;
				w.w[w.index(k,j,i)][1]=imag(aux)*factor;
			}       
	/**********************/
	fftw_execute(w.p3DB);
	/**********************/     
}






//Potential Evolver Version 1
void prop_potential(wavefunction &w, hamiltonian h, complex dt )
{
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				complex potoperator=exp(-I*dt*h.v[h.index(k,j,i)] );
				
				double aux0r=w.w[w.index(k,j,i)][0];
				double aux0i=w.w[w.index(k,j,i)][1];
				
				w.w[w.index(k,j,i)][0]=aux0r*real(potoperator)+aux0i*imag(potoperator);
				w.w[w.index(k,j,i)][1]=aux0i*real(potoperator)-aux0r*imag(potoperator);
			}
}



//Potential Evolver Version 2
void prop_potential_complex(wavefunction &w, hamiltonian h, double dt ) //For propagation in imaginary time
{
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				complex potoperator=complex(exp(-dt*h.v[h.index(k,j,i)] ), 0);
				complex aux=complex(w.w[w.index(k,j,i)][0], w.w[w.index(k,j,i)][1] );
				aux=aux*potoperator;
				
				w.w[w.index(k,j,i)][0]=real(aux);
				w.w[w.index(k,j,i)][1]=imag(aux);
			}
}




//Potential Evolver Version 3
void prop_potential_change(wavefunction &w, hamiltonian h, double dt ) //For propagation in real time
{
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				complex potoperator	=complex(cos(dt*h.v[h.index(k,j,i)] ), -sin(dt*h.v[h.index(k,j,i)] ) );
				complex aux			=complex(w.w[w.index(k,j,i)][0],w.w[w.index(k,j,i)][1])*potoperator;
				
				w.w[w.index(k,j,i)][0]= real(aux);
				w.w[w.index(k,j,i)][1]= imag(aux);
			}
}





//Potential Evolver Version 4
void prop_potential_length_gauge(wavefunction &w, hamiltonian h, double dt, double field_e1, double field_e2, double field_e3 ) //For propagation in real time
{
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)	
			{
				complex potoperator		= complex( cos( dt*(h.v[h.index(k,j,i)] 
															- h.q1*h.x1[i]*field_e1
															- h.q1*h.x2[j]*field_e2
															- h.q1*h.x3[k]*field_e3) ), 
												  -sin( dt*(h.v[h.index(k,j,i)] 
															- h.q1*h.x1[i]*field_e1
															- h.q1*h.x2[j]*field_e2
															- h.q1*h.x3[k]*field_e3) ) );
				
				complex aux				= potoperator*complex(w.w[w.index(k,j,i)][0], w.w[w.index(k,j,i)][1]);
				
				w.w[w.index(k,j,i)][0]	= real(aux);
				w.w[w.index(k,j,i)][1]	= imag(aux);
			}
}



//Absorber wavefunction 
void absorber(wavefunction w, double _frac_x1_left,double _frac_x2_left,double _frac_x3_left,double _frac_x1_right,double _frac_x2_right, double _frac_x3_right, double _exponent)
{    
	double frac_x1_right=_frac_x1_right;
	double frac_x2_right=_frac_x2_right;
	double frac_x3_right=_frac_x3_right;
	
	double frac_x1_left=_frac_x1_left;
	double frac_x2_left=_frac_x2_left;
	double frac_x3_left=_frac_x3_left;
	
	double exponent=_exponent; //try use 1./6.
	
	double mask_start_x1_right=w.x1[int(w.n1*(1.-frac_x1_right))];
	double mask_start_x2_right=w.x2[int(w.n2*(1.-frac_x2_right))];
	double mask_start_x3_right=w.x3[int(w.n3*(1.-frac_x3_right))];
	
	double mask_start_x1_left=w.x1[int(w.n1*frac_x1_left)+1];
	double mask_start_x2_left=w.x2[int(w.n2*frac_x2_left)+1];
	double mask_start_x3_left=w.x3[int(w.n3*frac_x3_left)+1];
	
	double argument_x1_right;
	double argument_x2_right;
	double argument_x3_right;
	
	double argument_x1_left;
	double argument_x2_left;	
	double argument_x3_left;
	
	double mask_x1_right;
	double mask_x2_right;
	double mask_x3_right;
	
	double mask_x1_left;
	double mask_x2_left;
	double mask_x3_left;
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				argument_x1_right=(pi/2.)*(w.x1[i]-mask_start_x1_right)/(w.x1[w.n1-1]-mask_start_x1_right+1.e-20);
				argument_x2_right=(pi/2.)*(w.x2[j]-mask_start_x2_right)/(w.x2[w.n2-1]-mask_start_x2_right+1.e-20);
				argument_x3_right=(pi/2.)*(w.x3[k]-mask_start_x3_right)/(w.x3[w.n3-1]-mask_start_x3_right+1.e-20);
				
				argument_x1_left=(pi/2.)*(w.x1[i]-mask_start_x1_left)/(w.x1[0]-mask_start_x1_left+1.e-20);
				argument_x2_left=(pi/2.)*(w.x2[j]-mask_start_x2_left)/(w.x2[0]-mask_start_x2_left+1.e-20);
				argument_x3_left=(pi/2.)*(w.x3[k]-mask_start_x3_left)/(w.x3[0]-mask_start_x3_left+1.e-20);
				
				mask_x1_right=pow(fabs(cos(argument_x1_right)),exponent);
				mask_x2_right=pow(fabs(cos(argument_x2_right)),exponent);
				mask_x3_right=pow(fabs(cos(argument_x3_right)),exponent);
				
				mask_x1_left=pow(fabs(cos(argument_x1_left)),exponent);
				mask_x2_left=pow(fabs(cos(argument_x2_left)),exponent);
				mask_x3_left=pow(fabs(cos(argument_x3_left)),exponent);
				
				if (i< int(w.n1*frac_x1_left))
				{		  
					w.w[w.index(k,j,i)][0]*=mask_x1_left;	
					w.w[w.index(k,j,i)][1]*=mask_x1_left;	
				}
				
				if (i> int(w.n1*(1.-frac_x1_right)))
				{		  
					w.w[w.index(k,j,i)][0]*=mask_x1_right;	
					w.w[w.index(k,j,i)][1]*=mask_x1_right;	
				}
				
				if (j< int(w.n2*frac_x2_left))
				{		  
					w.w[w.index(k,j,i)][0]*=mask_x2_left;	
					w.w[w.index(k,j,i)][1]*=mask_x2_left;	
				}
				
				if (j> int(w.n2*(1.-frac_x2_right)))
				{		  
					w.w[w.index(k,j,i)][0]*=mask_x2_left;	
					w.w[w.index(k,j,i)][1]*=mask_x2_left;	
				}
				
				if (k< int(w.n3*frac_x3_left))
				{		  
					w.w[w.index(k,j,i)][0]*=mask_x3_left;	
					w.w[w.index(k,j,i)][1]*=mask_x3_left;	
				}
				
				if (k> int(w.n3*(1.-frac_x3_right)))
				{		  
					w.w[w.index(k,j,i)][0]*=mask_x3_left;	
					w.w[w.index(k,j,i)][1]*=mask_x3_left;	
				}
			}//End the loop on kji
	
}



//Kinetic energy expectation value version 1
double kinetic_finite_diff(wavefunction w, hamiltonian h)
{
	complex kinetic2=0.;
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				complex phi=complex(w.w[w.index(k,j,i)][0],w.w[w.index(k,j,i)][1]);
				complex phi_ip=complex(0.,0.);
				complex phi_im=complex(0.,0.);
				complex phi_jp=complex(0.,0.); 
				complex phi_jm=complex(0.,0.);	     
				complex phi_kp=complex(0.,0.); 
				complex phi_km=complex(0.,0.);
				
				/******************************************************************/
				if (i-1>=0)
					phi_ip=complex(w.w[w.index(k,j,i-1)][0],w.w[w.index(k,j,i-1)][1]);		
				
				if(i+1<w.n1)
					phi_im=complex(w.w[w.index(k,j,i+1)][0],w.w[w.index(k,j,i+1)][1]);
				
				/******************************************************************/
				if (j-1>=0)
					phi_jp=complex(w.w[w.index(k,j-1,i)][0],w.w[w.index(k,j-1,i)][1]);
				
				if(j+1<w.n2)
					phi_jm=complex(w.w[w.index(k,j+1,i)][0],w.w[w.index(k,j+1,i)][1]);
				
				/******************************************************************/
				if (k-1>=0)
					phi_kp=complex(w.w[w.index(k-1,j,i)][0],w.w[w.index(k-1,j,i)][1]);		
				
				if(k+1<w.n3)
					phi_km=complex(w.w[w.index(k+1,j,i)][0],w.w[w.index(k+1,j,i)][1]);
				
				/******************************************************************/
				
				complex ksquare=
				-h.a1*(phi_ip-2.*phi+phi_im)/w.dx1/w.dx1
				-h.a2*(phi_jp-2.*phi+phi_jm)/w.dx2/w.dx2
				-h.a3*(phi_kp-2.*phi+phi_km)/w.dx3/w.dx3;
				
				kinetic2+=w.dx1*w.dx2*w.dx3*(conj(phi)*ksquare);
				
			}
	return real(kinetic2);
}


//Potential Energy Observable: Expectation potential energy value //
double potential_energy(wavefunction w, hamiltonian h)
{  
	double potE=0.;
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				potE+=w.dx1*w.dx2*w.dx3*h.v[h.index(k,j,i)]*
				(w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1]);
			}
	return potE;
}


//Potential Energy Observable: Expectation potential energy value into dynamic electric field //
double potential_energy_length_gauge1D(wavefunction w, hamiltonian h,double efield1)
{  
	double potE=0.;
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				potE+=w.dx1*w.dx2*w.dx3*(h.v[h.index(k,j,i)]+efield1*h.x1[i])*
				(w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1]);
			}
	return potE;
}


//Momentum distribution observable 
void momentum_distribution(wavefunction &w , wavefunction &momentum_w , double x0, double y0, double z0)
{
	
	double r;
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				if(w.n3==1 && w.dx3==1.)
					r=sqrt((w.x1[i]-x0)*(w.x1[i]-x0)+(w.x2[j]-y0)*(w.x2[j]-y0));
				else
					r=sqrt((w.x1[i]-x0)*(w.x1[i]-x0)+(w.x2[j]-y0)*(w.x2[j]-y0)+(w.x3[k]-z0)*(w.x3[k]-z0));	  
				double r0=10.;   //5.;
				double r1=20.;   //15.;
				
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/16.));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/16.));
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
				}
				
			}
    //momentum_w.normalize();
    fftw_execute(momentum_w.p3DF);
}



void Mask1D( wavefunction &w1, wavefunction &momentum_w, double _r0, double _r1)
{
	
	double r0=_r0;//10.;//5.;
	double r1=_r1;//20.;//15.;
	
	double r;
	for(int k=0;k<w1.n3;k++)
		for(int j=0;j<w1.n2;j++)
			for(int i=0;i<w1.n1;i++)
			{
				r=abs(w1.x1[i]);//sqrt( (w1.x1[i]-x0)*(w1.x1[i]-x0) );
				double sigma=abs(r1-r0)/5.;
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w1.w[w1.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w1.w[w1.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w1.w[w1.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w1.w[w1.index(k,j,i)][1];
				}
				
			}	
}



//***********Inverse mask function 1D *********//
void IMask1D( wavefunction &wfull, wavefunction &wmask,  double _r0, double _r1)
{
	
	double r0=_r0;//10.;//5.;
	double r1=_r1;//20.;//15.;
	
	double r;
	for(int k=0;k<wmask.n3;k++)
		for(int j=0;j<wmask.n2;j++)
			for(int i=0;i<wmask.n1;i++)
			{
				r=abs(wmask.x1[i]);//sqrt( (w1.x1[i]-x0)*(w1.x1[i]-x0) );
				double sigma=abs(r1-r0)/5.;
				if(r<=r0)
				{
					wmask.w[wmask.index(k,j,i)][0]=wfull.w[wfull.index(k,j,i)][0];
					wmask.w[wmask.index(k,j,i)][1]=wfull.w[wfull.index(k,j,i)][1];
				}
				if(r>r0 && r<=r1 )
				{
					wmask.w[wmask.index(k,j,i)][0]=wfull.w[wfull.index(k,j,i)][0]*exp(-(r-r0)*(r-r0)/sigma/sigma);
					wmask.w[wmask.index(k,j,i)][1]=wfull.w[wfull.index(k,j,i)][1]*exp(-(r-r0)*(r-r0)/sigma/sigma);
				}
				if(r>r1)
				{
					wmask.w[wmask.index(k,j,i)][0]=0.;
					wmask.w[wmask.index(k,j,i)][1]=0.;
				}
				
			}	
}



//***********Inverse mask function 1D *********//
void IMask2D( wavefunction &wfull, wavefunction &wmask,  double _r0, double _r1)
{
	
	double r0=_r0;//10.;//5.;
	double r1=_r1;//20.;//15.;
	double sigma=abs(r1-r0)/5.;	
	
	double r;
	
	
	for(int k=0;k<wmask.n3;k++)
		for(int j=0;j<wmask.n2;j++)
			for(int i=0;i<wmask.n1;i++)
			{
				r=sqrt(wmask.x1[i]*wmask.x1[i] + wmask.x2[j]*wmask.x2[j]);//sqrt( (w1.x1[i]-x0)*(w1.x1[i]-x0) );

				if(r<=r0)
				{
					wmask.w[wmask.index(k,j,i)][0]=wfull.w[wfull.index(k,j,i)][0];
					wmask.w[wmask.index(k,j,i)][1]=wfull.w[wfull.index(k,j,i)][1];
				}
				if(r>r0 && r<=r1 )
				{
					wmask.w[wmask.index(k,j,i)][0]=wfull.w[wfull.index(k,j,i)][0]*exp(-(r-r0)*(r-r0)/sigma/sigma);
					wmask.w[wmask.index(k,j,i)][1]=wfull.w[wfull.index(k,j,i)][1]*exp(-(r-r0)*(r-r0)/sigma/sigma);
				}
				if(r>r1)
				{
					wmask.w[wmask.index(k,j,i)][0]=0.;
					wmask.w[wmask.index(k,j,i)][1]=0.;
				}
				
			}	
}



//Mask functions //
void Mask_Shift1D( wavefunction &w1, wavefunction &momentum_w, double _r0, double _r1, double R0 )
{
	
	double r0=_r0;//10.;//5.;
	double r1=_r1;//20.;//15.;
	
	double r;
	for(int k=0;k<w1.n3;k++)
		for(int j=0;j<w1.n2;j++)
			for(int i=0;i<w1.n1;i++)
			{
				r=abs(w1.x1[i]-R0);//sqrt( (w1.x1[i]-x0)*(w1.x1[i]-x0) );
				double sigma=abs(r1-r0)/5.;
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w1.w[w1.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w1.w[w1.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w1.w[w1.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w1.w[w1.index(k,j,i)][1];
				}
				
			}	
}


// Mask functions //
void Mask2D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1)
{
	
	double r;
	for(int k=0;k<w.n3;k++)	
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				r=sqrt( w.x1[i]*w.x1[i]+w.x2[j]*w.x2[j] );
				
				double r0=_r0;//5.;
				double r1=_r1;//15.;
				double sigma=(r1-r0)/5.;
				//double r0=20.;
				//double r1=30.;
				
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));;
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
				}
				
			}
	
}



//Anti Mask function 
void AntiMask2D( wavefunction &w1 , double _r0, double _r1, double _sigma0 )
{
	
	double r0=_r0;//10.;//5.;
	double r1=_r1;//20.;//15.;
	
	double r;
	for(int k=0;k<w1.n3;k++)
		for(int j=0;j<w1.n2;j++)
			for(int i=0;i<w1.n1;i++)
			{
				r=sqrt( w1.x1[i]*w1.x1[i] + w1.x2[j]*w1.x2[j]  );
				if(r>r0 && r<=r1 )
				{
					w1.w[w1.index(k,j,i)][0]*=(exp(-(r-r0)*(r-r0)/_sigma0));
					w1.w[w1.index(k,j,i)][1]*=(exp(-(r-r0)*(r-r0)/_sigma0));;
				}
				if(r>r1)
				{
					w1.w[w1.index(k,j,i)][0]=0.;
					w1.w[w1.index(k,j,i)][1]=0.;
				}
				
			}
	
}



//Momentum distribution observable 1D
void momentum_distribution1D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1)
{
	
	double r;
	set_wavefunction_zeros(momentum_w);
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				r=sqrt( (w.x1[i])*(w.x1[i]) );
				double r0=_r0;//10.;//5.;
				double r1=_r1;//30.;//15.;
				double sigma=abs(r1-r0)/5;
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));     //50.));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));     //50.));
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
				}				
			}
	//momentum_w.wzero(pow(10.,-100.), pow(10.,-100.));
	fftw_execute(momentum_w.p3DF);	
}

//Momentum distribution observable 2D
void momentum_distribution2D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1 )
{	
	double r;
	for(int k=0;k<w.n3;k++)	
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				r=sqrt( w.x1[i]*w.x1[i]+w.x2[j]*w.x2[j] );
				
				double r0=_r0;//5.;
				double r1=_r1;//15.;
				double sigma=(r1-r0)/5.;
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/sigma/sigma));;
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
				}
				
			}
	fftw_execute(momentum_w.p3DF);	
}



//Mask function 1D
void mask_function1D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1, double _sigma )
{
	
	double r;
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				r=sqrt( (w.x1[i])*(w.x1[i]) );
				
				double r0=_r0;//10.;//5.;
				double r1=_r1;//30.;//15.;
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/_sigma));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/_sigma));;
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
				}				
			}
}



//Kinetic Energy obsevable: Kinetic Energy Expectation value //
double q_expected_kinetic(wavefunction &w, hamiltonian &h)
{
	w.expected_kin=0.;
	w.qnorm=0.;
	double factor=1./w.n1/w.n2/w.n3;
	
	/**********************/
	fftw_execute(w.p3DF);
	/**********************/
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				
				double kinetic=
				w.q1[i]*w.q1[i]*h.a1+
				w.q2[j]*w.q2[j]*h.a2+
				w.q3[k]*w.q3[k]*h.a3;
				
				w.expected_kin+=w.dq1*w.dq2*w.dq3*kinetic*
				(w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1])
				*w.kfactor1*w.kfactor1*w.kfactor2*w.kfactor2*w.kfactor3*w.kfactor3;
				
				w.w[w.index(k,j,i)][0]*=factor;
				w.w[w.index(k,j,i)][1]*=factor;
				
			}    
	/**********************/
	fftw_execute(w.p3DB);
	/**********************/
	
	return w.expected_kin;
}



//Second Version of Kinetic Energy obsevable: Kinetic Energy Expectation value //
double q_expected_kinetic_copy(wavefunction &w, hamiltonian &h)
{  
	grid *g1=new grid[1];  
	g1[0].set_grid(w.n1,w.n2,w.n3,w.dx1,w.dx2,w.dx3);
	
	wavefunction *wcopy=new wavefunction[1];
	wcopy[0].put_on_grid(g1[0]);
	
	placeWF(wcopy[0], w);
	
	wcopy[0].expected_kin=0.;
	wcopy[0].expected_kin=q_expected_kinetic(wcopy[0], h);
	
	return wcopy[0].expected_kin;
	
	delete[] g1;
	delete[] wcopy;
}



//This function calculate the complex number <w2|w1> //
complex projection(wavefunction &w2, wavefunction &w1)
{
	
	complex projection =complex(0.,0.);
	for(int k=0;k<w1.n3;k++)
		for(int j=0;j<w1.n2;j++)
			for(int i=0;i<w1.n1;i++)
			{
				complex wave1=complex(w1.w[w1.index(k,j,i)][0],w1.w[w1.index(k,j,i)][1]);
				complex wave2=complex(w2.w[w2.index(k,j,i)][0],w2.w[w1.index(k,j,i)][1]);
				
				projection+=w1.dx1*w1.dx2*w1.dx3*conj(wave2)*wave1;
			}
	
	return projection;
}



//Scattering wave projector 1D //
complex scattering_wave_projector1D(wavefunction &w, hamiltonian &h, wavefunction &wscattering, double kx )
{
	
	Continuum_WF( wscattering, h, kx);		
	complex proj = projection( wscattering, w );
	
	
	return proj;
	
}


//Set of Scattering wave projector 1D //
void set_scattering_wave_projector1D(wavefunction &w, hamiltonian &h, complex *proj, int Nmom, double dk )
{
	//Parallel OpenMP version 	
#pragma omp parallel shared(proj,dk,Nmom)	
		{		
			grid *g=new grid[1];
			g[0].set_grid(w.n1, w.n2, w.n3, w.dx1, w.dx2, w.dx3);
			
			
			wavefunction *wscattering=new wavefunction[1];
			wscattering[0].put_on_grid(g[0]);
			
			int imom;
			double kmin =-(Nmom-1)/2*dk;
			double kx	=0;
			
#pragma omp for private(imom,kx)	
			for (imom=0; imom<Nmom; imom++) 
			{
				kx			= kmin + dk*imom;
				proj[imom]	= scattering_wave_projector1D(w, h, wscattering[0], kx );
				
			}//End Projection
			delete[] wscattering;
			delete[] g;						
		}//End pragma		*/
	

}


//Excited population 
double bound_population1D(wavefunction &w, wavefunction *warray, int Nstates, int index_initial_state)
{
	double temp0=0.;
	complex temp1;
	for (int i=0; i<Nstates; i++){
		temp1 = projection(warray[i], w);
		if(i!=index_initial_state)
			temp0+=abs(temp1)*abs(temp1);
	}

	return temp0;
}


//Ionization Population 1D By Scattering Wave projection 
double ionization_scatteringW1D(complex *proj, int Nmom, double dk)
{
	double temp0=0.;
	for(int imom=0;imom<Nmom;imom++)
		temp0+=abs(proj[imom])*abs(proj[imom])*dk;
	
	return temp0;
}



//Projection <w2|w1> 
complex projection_on_range(wavefunction &w2, wavefunction &w1,double *x,double *y,double *z)
{
    //Projection on right and left direction
	//if sw0 =0, projection on left
	//if sw0 !=0 projection on right
	int index_x0 = floor((x[0] - w1.x1[0])/w1.dx1);
	int index_x1 = floor((x[1] - w1.x1[0])/w1.dx1)+1;	
	
	if (index_x0<0) {
		cout << "\n\n//*******TAKE CARE:\nThe variable x[0] must be bigger than grid x1[0] point\n\n ";
		exit(1);
	}
	
	//Checking y dimension
	int index_y0=0;
	int index_y1=1; 	
	if(w1.n2>1){
		index_y0	= floor((y[0] - w1.x2[0])/w1.dx2);
		index_y1	= floor((y[1] - w1.x2[0])/w1.dx2)+1;		
	}
	
	//Checking z dimension
	int index_z0 = 0;
	int index_z1 = 1;		
	
	if (w1.n3>1){
		index_z0 = floor((z[0] - w1.x3[0])/w1.dx3);
		index_z1 = floor((z[1] - w1.x3[0])/w1.dx3)+1;		
	}

	
	
	complex projection =complex(0.,0.);

	for(int k=index_z0;k<index_z1;k++)
		for(int j=index_y0;j<index_y1;j++)
			for(int i=index_x0;i<index_x1;i++)
			{
				complex wave1=complex(w1.w[w1.index(k,j,i)][0],w1.w[w1.index(k,j,i)][1]);
				complex wave2=complex(w2.w[w2.index(k,j,i)][0],w2.w[w1.index(k,j,i)][1]);
				
				projection+=w1.dx1*w1.dx2*w1.dx3*conj(wave2)*wave1;
			}

	return projection;
}

//Position Population On Range //
double position_population(wavefunction &w, double *x, double *y, double *z)
{
	complex cpop = projection_on_range(w,w,x,y,z);
	
	return abs(cpop);
	
}



//Momentum Population On Range (x[0] x[1]), (y[0] y[1]), (z[0] z[1]) //
/*double momentum_population1D(wavefunction &w, double qxa, double qxb)
{	
	
	if(qxa>qxb){
		cout << "\n\n****Please check the integation limits qxa and qxb \n\n";
		exit(1);
	}
	
	w.expected_kin=0.;
	w.qnorm=0.;
	double factor	=1./w.n1;
	int Nqa			= 0;
	int Nqb			= 0;
	

	//fftw_execute(w.p3DF);


	if (qxa<0. && qxb>0.){
	
		Nqa = floor()+;
		
		for(int i=0;i<w.n1;i++)
			{	  
				double kinetic=
				w.q1[i]*w.q1[i]*h.a1+
				w.q2[j]*w.q2[j]*h.a2+
				w.q3[k]*w.q3[k]*h.a3;
				
				complex kinoperator=exp(-I*dt*kinetic );
				
				double aux0r=w.w[w.index(k,j,i)][0];
				double aux0i=w.w[w.index(k,j,i)][1];
				
				w.w[w.index(k,j,i)][0]=( aux0r*real(kinoperator)+aux0i*imag(kinoperator) )*factor;
				w.w[w.index(k,j,i)][1]=( aux0i*real(kinoperator)-aux0r*imag(kinoperator) )*factor;
			}  
		
	}

	fftw_execute(w.p3DB);

	
	return abs(cpop);
	
}*/



//*****Remove w2 from w1
void project_out(wavefunction &w1, wavefunction &w2)
{
	
	complex p2=projection(w2,w1);
	
	for(int k=0;k<w1.n3;k++)
		for(int j=0;j<w1.n2;j++)
			for(int i=0;i<w1.n1;i++)
			{
				complex wave1=complex(w1.w[w1.index(k,j,i)][0],w1.w[w1.index(k,j,i)][1]);
				complex wave2=complex(w2.w[w2.index(k,j,i)][0],w2.w[w1.index(k,j,i)][1]);
				
				complex w=wave1-p2*wave2;
				w1.w[w1.index(k,j,i)][0]=real(w);
				w1.w[w1.index(k,j,i)][1]=imag(w);
			}
}



//Position Space Snapshot Probability from the wavefunction//
void snapshot(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3)
{
	
	for(int k=0;k<w1.n3/skiper3;k++)
		for(int j=0;j<w1.n2/skiper2;j++)
			for(int i=0;i<w1.n1/skiper1;i++)
			{
				
				double norm= w1.dx1*w1.dx2*w1.dx3*(
												   w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*
												   w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]+
												   w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*
												   w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]);
				fprintf(file,"%12.16e\n",norm);
			}
	fflush(file);
}


//Position Space Snapshot the wavefunction//
void complex_snapshot_WF(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3)
{
	
	for(int k=0;k<w1.n3/skiper3;k++)
		for(int j=0;j<w1.n2/skiper2;j++)
			for(int i=0;i<w1.n1/skiper1;i++)
				fprintf(file,"%e %e\n", w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0],
										w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]);

}


//Momentum Space Snapshot wavefunction//
void Qsnapshot(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3)
{
	
	for(int k=0;k<w1.n3/skiper3;k++)
		for(int j=0;j<w1.n2/skiper2;j++)
			for(int i=0;i<w1.n1/skiper1;i++)
			{
				double qnorm=w1.dq1*w1.dq2*w1.dq3*
				(w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0] +
				 
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1])
				*w1.kfactor1*w1.kfactor1*w1.kfactor2*w1.kfactor2*w1.kfactor3*w1.kfactor3;
				fprintf(file,"%e\n",qnorm);
			}
}




//Save MOMENTUM DISTRIBUTION 
void Qsnapshot_WF(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3)
{
	
	double factor=1./ w1.n1/ w1.n2/ w1.n3;
	double konst=w1.kfactor1*w1.kfactor2*w1.kfactor3;
	w1.Direct_Fourier_Transform();	

	for(int k=0;k<w1.n3/skiper3;k++)
		for(int j=0;j<w1.n2/skiper2;j++)
			for(int i=0;i<w1.n1/skiper1;i++)
			{
				double qnorm=w1.dq1*w1.dq2*w1.dq3*
				(w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0] +
				 
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1])*konst*konst;
				
				fprintf(file,"%e\n",qnorm);
				
				
			}
	fflush(file);
	cproduct( w1, complex(factor,0.));
	w1.Back_Fourier_Transform();
}


//Complex QSnapshot wavefunction
void complex_Qsnapshot_WF(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3)
{
	double factor=1./ w1.n1/ w1.n2/ w1.n3;
	double konst=w1.kfactor1*w1.kfactor2*w1.kfactor3;
	w1.Direct_Fourier_Transform();	
	
	for(int k=0;k<w1.n3/skiper3;k++)
		for(int j=0;j<w1.n2/skiper2;j++)
			for(int i=0;i<w1.n1/skiper1;i++)
				fprintf(file,"%e %e\n", konst*w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0],
										konst*w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]);

	fflush(file);
	cproduct( w1, complex(factor,0.));
	w1.Back_Fourier_Transform();	
	
}


//Complex amplitude transition by projection on scattering wave
void esay_complex_snapshot_scatteringW1D(FILE *file, wavefunction &w, hamiltonian &h, int Nmom, double dk)
{
	complex *proj;
	proj		= (complex*)malloc(Nmom*sizeof(complex));
	
	complex_snapshot_scatteringW1D(file, w, h, proj, Nmom, dk);
	
	
	free(proj);
}

//Complex amplitude transition by projection on scattering wave
void complex_snapshot_scatteringW1D(FILE *file, wavefunction &w, hamiltonian &h, complex *proj, int Nmom, double dk)
{
	
	//Set Scattering Wave Projector 
	set_scattering_wave_projector1D(w, h, proj, Nmom, dk );
	
	//Saving Data
	for(int imom=0;imom<Nmom;imom++)
		fprintf(file,"%e %e\n",real(proj[imom]),imag(proj[imom]));
	fflush(file);
	
}



//Square Amplitude transition by projection on scattering wave
void easy_snapshot_scatteringW1D(FILE *file, wavefunction &w, hamiltonian &h,  int Nmom, double dk)
{
	
	complex *proj;
	proj		= (complex*)malloc(Nmom*sizeof(complex));
	snapshot_scatteringW1D( file, w, h, proj, Nmom, dk);
	free(proj);
	
}


//Square Amplitude transition by projection on scattering wave
void snapshot_scatteringW1D(FILE *file, wavefunction &w, hamiltonian &h, complex *proj, int Nmom, double dk)
{
	//Set Scattering Wave Projector 
	set_scattering_wave_projector1D(w, h, proj, Nmom, dk );
	
	for(int imom=0;imom<Nmom;imom++)
		fprintf(file,"%e\n",abs(proj[imom])*abs(proj[imom])*dk);
	fflush(file);	
}


//Square Amplitude transition bound to bound states
void snapshot_bound_to_bound_transition(FILE *file, wavefunction &w, wavefunction *warray, int Nstates)
{
	
	complex temp0;
	for(int i=0;i<Nstates;i++){
		temp0 = projection(warray[i],w);		
		fprintf(file,"%e\n",abs(temp0)*abs(temp0));
	}
	fflush(file);
}


//Binary wave save in position space
void biwrite(FILE *file, wavefunction &w)  
{
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++){
				
				double wreal= w.w[w.index(k,j,i)][0];
				double wimag= w.w[w.index(k,j,i)][1];
				
				fwrite (&wreal , 1 , sizeof(wreal) , file );
				fwrite (&wimag , 1 , sizeof(wimag) , file );
			}
	fflush(file);
}



//Binary reader in position space
void biread(FILE *file, wavefunction &w)   
{
	double wreal;
	double wimag;
	
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++){
				
				fread(&wreal , 1 , sizeof(wreal) , file );
				fread(&wimag , 1 , sizeof(wimag) , file );
				
				w.w[w.index(k,j,i)][0]=wreal;
				w.w[w.index(k,j,i)][1]=wimag;
			}    
	fflush(file);	
}




//Saving axis, position and momentum //
void XQsnapshot(FILE *file, grid &g1, int skiper1, int skiper2, int skiper3)
{
	if (g1.n1>1 )
	{
		for(int i=0;i<g1.n1/skiper1;i++) 	    
	        fprintf(file,"%e %e\n", g1.x1[g1.index1_1d(i*skiper1)], g1.q1[g1.index1_1d(i*skiper1)]);
	}
	if (g1.n2>1 )
	{
		for(int j=0; j<g1.n2/skiper2; j++) 	    
	        fprintf(file,"%e %e\n", g1.x2[g1.index2_1d(j*skiper2)], g1.q2[g1.index2_1d(j*skiper2)]);
	}
	if (g1.n3>1 )
	{
		for(int k=0; k<g1.n3/skiper3; k++) 	    
	        fprintf(file,"%e %e\n", g1.x3[g1.index3_1d(k*skiper3)], g1.q3[g1.index3_1d(k*skiper3)]);
	}
}




//Reading field from .txt format 
double wlecture(FILE *file, wavefunction &w)
{
	double temp0;
	double temp1;
	double eng0;
	double eng;
	for(int i=0;i<w.n1*w.n2*w.n3+1;i++)
	{
		if (i==0)
		{
			fscanf(file,"%lf %lf\n",&temp0,&temp1);
			eng0=temp0;
			eng=temp1;
		}
		else
		{
			fscanf(file,"%lf %lf\n",&temp0,&temp1);
			w.w[i-1][0]=temp0;
			w.w[i-1][1]=temp1;
		}
	}
	return eng0;
}




//Definite integral by trapezoidal method
double trapz_1D(vector <double> &x, vector<double> &y, double a, double b)
{	/*****************============================
	//Definite integral by trapezoidal method
	//For the function y(x) from a to b
	//**************===============================  */
	
	int nx        = x.size();
	double dx     = x[1]-x[0];	
	double answer = 0.;
	
	int Na = floor( abs( ( a - x[0] )/dx  )  );
	int Nb = floor( abs( ( b - x[0] )/dx  )  ) + 1;
	
	if ((b-a)>0.) {
		for (int i=Na; i<Nb-1; i++)
			answer+= dx*( y[i] + y[i+1] )/2.;
	}
	else {
		dx=-dx;
		for (int i=Nb-1; i<Na; i++) 
			answer+= dx*( y[i] + y[i+1] )/2.;
	}

	return answer;	
}




//This routine make the integral of the ...
complex def_integrate( wavefunction &w1, double a1, double b1, double a2, double b2, double a3, double b3)
{
	//This routine make the integral of the
	//wavefunction from a0 to b0  x-axis
	//wavefunction from a1 to b1  y-axis
	//wavefunction from a2 to b2  z-axis
	//Introduce a1=b1 and a2=b2 to make the integration though x-direction
	//Introduce a0=b0 and a2=b2 to make the integration though y-direction
	//Introduce a0=b0 and a1=b1 to make the integration though z-direction	
	//...
	
	complex int0=complex(0.,0.);
	double s1=1.0;
	double s2=1.0;
	double s3=1.0;	

	
	if (a1==b1) s1=0.;
	if (a2==b2) s2=0.;
	if (a3==b3) s3=0.;	

	
	int ka1 = floor(s1*abs(a1 - w1.x1[0])/w1.dx1  ) ;
	int ka2 = floor(s2*abs(a2 - w1.x2[0])/w1.dx2  ) ;
	int ka3 = floor(s3*abs(a3 - w1.x3[0])/w1.dx3  ) ;
	
	
	int kb1 = floor(s1*abs(b1 - w1.x1[0])/w1.dx1  ) +1;
	int kb2 = floor(s2*abs(b2 - w1.x2[0])/w1.dx2  ) +1;
	int kb3 = floor(s3*abs(b3 - w1.x3[0])/w1.dx3  ) +1;
	
	if (b1-a1>0) {

		for(int k=ka3; k<kb3; k++)
			for(int j=ka2; j<kb2; j++)
				for(int i=ka1; i<kb1; i++)				
					int0+= complex( w1.w[w1.index(k,j,i)][0], w1.w[w1.index(k,j,i)][1] )*w1.dx1*w1.dx2*w1.dx3;
	}
	else {
		for(int k=ka3; k<kb3; k++)
			for(int j=ka2; j<kb2; j++)
				for(int i=kb1-1; i<ka1; i++)		{
					int0+= complex( w1.w[i][0], w1.w[i][1] )*w1.dx1*w1.dx2*w1.dx3;//complex( w1.w[ka1-i][0], w1.w[ka1-i][1] )*w1.dx1*w1.dx2*w1.dx3;
					
					//if (i%1==0) 
					//	cout << "\ni:  "<<i;//ka1-i;
					

				}
		
	}


	return int0;
}





//Linear combination: w3=b[0]*w1+b[1]*w2; b[0] and b[1] are reals
void linear_combination(wavefunction &w1, wavefunction &w2, wavefunction &w3, vector<double> &b )
{
    if (w1.n1*w1.n2*w1.n3!=w2.n1*w2.n2*w2.n3 || w1.n1*w1.n2*w1.n3!=w3.n1*w3.n2*w3.n3){
        cout<<"\nError in the size of the wavefunctions: w1, w2 and w3. "<<endl;
        cout<<"The size of the wavefunctions must be equal. "<<endl;        
	}          
	
	
	
	set_wavefunction_zeros(w3);
	
	for (int i=0;i<w1.n1*w1.n2*w1.n3;i++){
		w3.w[i][0] = b[0]*w1.w[i][0] + b[1]*w2.w[i][0];
		w3.w[i][1] = b[0]*w1.w[i][1] + b[1]*w2.w[i][1];
	}
	
	w3.normalize_q();
}



//Linear combination: w3=b[0]*w1+b[1]*w2; b[0] and b[1] are complex
void complex_linear_combination(wavefunction &w1, wavefunction &w2, wavefunction &w3, vector<complex> &b )
{
    if (w1.n1*w1.n2*w1.n3!=w2.n1*w2.n2*w2.n3 || w1.n1*w1.n2*w1.n3!=w3.n1*w3.n2*w3.n3){
        cout<<"Error in the size of the wavefunctions: w1, w2 and w3. "<<endl;
        cout<<"The size of the wavefunctions must be equal. "<<endl;        
	} 
	
	set_wavefunction_zeros(w3);
	
	for (int i=0;i<w1.n1*w1.n2*w1.n3;i++){
		complex aux1=b[0]*w1.w[i][0] + b[1]*w2.w[i][0];
		complex aux2=b[0]*w1.w[i][1] + b[1]*w2.w[i][1];
		
		w3.w[i][0] = real(aux1)+real(aux2);
		w3.w[i][1] = imag(aux1)+imag(aux2);
	}
	
	//w3.normalize_q();
}




//1D ordering q momentum grid//
void q_order(grid &g)
{
	vector<double> q1aux;
	vector<double> q2aux;
	vector<double> q3aux;
	
	q1aux.resize(g.n1,0);
	q2aux.resize(g.n2,0);
	q3aux.resize(g.n3,0);
	
	if(g.n1!=1){
		for (int i=0; i<g.n1; i++){
			if (i<g.n1/2)
                q1aux[i] = g.q1[i+g.n1/2];            
			else
                q1aux[i] = g.q1[i-g.n1/2];            
		}
	}
	
	if(g.n2!=1){
		for (int j=0; j<g.n2; j++){
			if (j<g.n2/2)
                q2aux[j] = g.q2[j+g.n2/2];            
			else
                q2aux[j] = g.q2[j-g.n2/2];            
		}
	}
	if(g.n3!=1){
		for (int k=0; k<g.n3; k++){
			if (k<g.n3/2)
                q3aux[k] = g.q3[k+g.n3/2];            
			else
                q3aux[k] = g.q3[k-g.n3/2];            
		}
	}
	
	for (int i=0; i<g.n1; i++)
		g.q1[i]   = q1aux[i];
	for (int j=0; j<g.n2; j++)
		g.q2[j]   = q2aux[j];
	for (int k=0; k<g.n3; k++)
		g.q3[k]   = q3aux[k];
}




//1D ordering wavefunction //
void w_order1D(wavefunction &w, int sw0)
{    
	//sw0=1, order the phase and the fft (w) of the wave function w
	vector<double> wauxR;
	vector<double> wauxI;
	wauxR.resize(w.n1,0);
	wauxI.resize(w.n1,0);         
	
	if(sw0==1){
		for(int i=0; i<w.n1; i++){
            if(i<w.n1/2){
                wauxR[i] = w.w[i+w.n1/2][0];
                wauxI[i] = w.w[i+w.n1/2][1];
			}
            else{
                wauxR[i] = w.w[i-w.n1/2][0];
                wauxI[i] = w.w[i-w.n1/2][1];
			}
		}
		
		for(int i=0; i<w.n1; i++){
       		w.w[i][0] = wauxR[i];
       		w.w[i][1] = wauxI[i];
		}
		
		for(int i=0; i<w.n1; i++){
            if (i<w.n1/2)
                wauxR[i] = w.Qphase[i+w.n1/2];             
            else
                wauxR[i] = w.Qphase[i-w.n1/2];             
		}
		for(int i=0; i<w.n1; i++)
			w.Qphase[i]=wauxR[i];
	}
	else{
		for(int i=0; i<w.n1; i++){
            if(i<w.n1/2){
                wauxR[i] = w.w[i+w.n1/2][0];
                wauxI[i] = w.w[i+w.n1/2][1];
			}
            else{
                wauxR[i] = w.w[i-w.n1/2][0];
                wauxI[i] = w.w[i-w.n1/2][1];
			}
		}
		for(int i=0; i<w.n1; i++){
       		w.w[i][0] = wauxR[i] ;
       		w.w[i][1] = wauxI[i];
		}       
	}
}


//2D ordering wavefunction //
void w_order2D(wavefunction &w, int sw0)
{    
	vector<double> wauxR;
	vector<double> wauxI;
	wauxR.resize(w.n1*w.n2,0);
	wauxI.resize(w.n1*w.n2,0);
	
	for(int j=0; j<w.n2; j++)
        for(int i=0; i<w.n1; i++){
			if(i<w.n1/2 && j<w.n2_half()){
				wauxR[w.index12_2d(j,i)]   = w.w[w.index12_2d(j,i)+w.n1*w.n2/2][0];
				wauxI[w.index12_2d(j,i)]   = w.w[w.index12_2d(j,i)+w.n1*w.n2/2][1];
			}
            else{
				wauxR[w.index12_2d(j,i)]   = w.w[w.index12_2d(j,i)-w.n1*w.n2/2][0];
				wauxI[w.index12_2d(j,i)]   = w.w[w.index12_2d(j,i)-w.n1*w.n2/2][1];
			}
		}
	
	for (int i=0; i<w.n1*w.n2; i++){
		w.w[i][0] = wauxR[i];
		w.w[i][1] = wauxI[i];
	}
	
	if(sw0==1){
		for(int j=0; j<w.n2; j++)
			for(int i=0; i<w.n1; i++){
                if(i<w.n1/2 && j<w.n2_half()){
					wauxR[w.index12_2d(j,i)] = w.w[w.index12_2d(j,i)+w.n1*w.n2/2][0];
					wauxI[w.index12_2d(j,i)] = w.w[w.index12_2d(j,i)+w.n1*w.n2/2][1];
				}
                else{
					wauxR[w.index12_2d(j,i)] = w.w[w.index12_2d(j,i)-w.n1*w.n2/2][0];
					wauxI[w.index12_2d(j,i)] = w.w[w.index12_2d(j,i)-w.n1*w.n2/2][1];
				}
			}
		
		for(int i=0; i<w.n1; i++){
       		w.w[i][0] = wauxR[i];
       		w.w[i][1] = wauxI[i];
		}
	}
}


//1D disordering grid//
void q_disorder(grid &g)
{
	vector<double> q1aux;
	vector<double> q2aux;
	vector<double> q3aux;
	
	q1aux.resize(g.n1,0);
	q2aux.resize(g.n2,0);
	q3aux.resize(g.n3,0);
	
	if(g.n1!=1){
		for (int i=0; i<g.n1; i++){
			if (i<g.n1/2)
                q1aux[i+g.n1/2] = g.q1[i];            
			else
                q1aux[i-g.n1/2] = g.q1[i];            
		}
	}
	
	if(g.n2!=1){
		for (int j=0; j<g.n2; j++){
			if (j<g.n2/2)
                q2aux[j+g.n2/2] = g.q2[j];            
			else
                q2aux[j-g.n2/2] = g.q2[j];            
		}
	}
	
	if(g.n3!=1){
		for (int k=0; k<g.n3; k++){
			if (k<g.n3/2)
                q3aux[k+g.n3/2] = g.q3[k];            
			else
                q3aux[k-g.n3/2] = g.q3[k];            
		}
	}
	
	for (int i=0; i<g.n1; i++)
		g.q1[i]   = q1aux[i];
	for (int j=0; j<g.n2; j++)
		g.q2[j]   = q2aux[j];
	for (int k=0; k<g.n3; k++)
		g.q3[k]   = q3aux[k];
}


//1D disordering wavefunction //
void w_disorder1D(wavefunction &w, int sw0)
{    
	vector<double> wauxR;
	vector<double> wauxI;
	wauxR.resize(w.n1,0);
	wauxI.resize(w.n1,0);
	
	for(int i=0; i<w.n1; i++){
		if(i<w.n1/2){
			wauxR[i+w.n1/2]   = w.w[i][0];
			wauxI[i+w.n1/2]   = w.w[i][1];
		}
        else{
			wauxR[i-w.n1/2]   = w.w[i][0];
			wauxI[i-w.n1/2]   = w.w[i][1];
		}
	}
	
	for (int i=0; i<w.n1; i++){
		w.w[i][0] = wauxR[i];
		w.w[i][1] = wauxI[i];
	}
	
	if(sw0==1){
		for(int i=0; i<w.n1; i++){
            if (i<w.n1/2){
                wauxR[i+w.n1/2] = w.w[i][0];
                wauxI[i+w.n1/2] = w.w[i][1];
			}
            else{
                wauxR[i-w.n1/2] = w.w[i][0];
                wauxI[i-w.n1/2] = w.w[i][1];
			}
		}
		
		for(int i=0; i<w.n1; i++){
       		w.w[i][0] = wauxR[i] ;
       		w.w[i][1] = wauxI[i];
		}
	}
}


//2D disordering wavefunction //
void w_disorder2D(wavefunction &w, int sw0)
{    
	vector<double> wauxR;
	vector<double> wauxI;
	wauxR.resize(w.n1*w.n2,0);
	wauxI.resize(w.n1*w.n2,0);
	
	for(int j=0; j<w.n2; j++)
        for(int i=0; i<w.n1; i++){
			if(i<w.n1/2 && j<w.n2_half()){
				wauxR[w.index12_2d(j,i)+w.n1*w.n2/2]   = w.w[w.index12_2d(j,i)][0];
				wauxI[w.index12_2d(j,i)+w.n1*w.n2/2]   = w.w[w.index12_2d(j,i)][1];
			}
            else{
				wauxR[w.index12_2d(j,i)-w.n1*w.n2/2]   = w.w[w.index12_2d(j,i)][0];
				wauxI[w.index12_2d(j,i)-w.n1*w.n2/2]   = w.w[w.index12_2d(j,i)][1];
			}
		}
	
	for (int i=0; i<w.n1*w.n2; i++){
		w.w[i][0] = wauxR[i];
		w.w[i][1] = wauxI[i];
	}
	
	if(sw0==1){
		for(int j=0; j<w.n2; j++)
			for(int i=0; i<w.n1; i++){
                if(i<w.n1/2 && j<w.n2_half()){
					wauxR[w.index12_2d(j,i)+w.n1*w.n2/2] = w.w[w.index12_2d(j,i)][0];
					wauxI[w.index12_2d(j,i)+w.n1*w.n2/2] = w.w[w.index12_2d(j,i)][1];
				}
                else{
					wauxR[w.index12_2d(j,i)-w.n1*w.n2/2] = w.w[w.index12_2d(j,i)][0];
					wauxI[w.index12_2d(j,i)-w.n1*w.n2/2] = w.w[w.index12_2d(j,i)][1];
				}
			}
		
		for(int i=0; i<w.n1; i++){
       		w.w[i][0] = wauxR[i];
       		w.w[i][1] = wauxI[i];
		}
	}
}




//Maxima value of the w
double wmax_x(wavefunction &w)
{
   	vector<double> w1;
	w1.resize(w.n1*w.n2*w.n3,0.0);
	for(int i=0;i<w.n1*w.n2*w.n3;i++)
		w1[i]=(w.w[i][0]*w.w[i][0]+w.w[i][1]*w.w[i][1])*w.dx1*w.dx2*w.dx3;
   	return qmajor(w1);
}





double wmax_q(wavefunction &w)
{
   	vector<double> w1;
	w1.resize(w.n1*w.n2*w.n3,0.0);
	for(int i=0;i<w.n1*w.n2*w.n3;i++)
		w1[i]=(w.w[i][0]*w.w[i][0]+w.w[i][1]*w.w[i][1])*w.dq1*w.dq2*w.dq3;
   	return qmajor(w1);
}






/******************************************************
 Numerov method plane wave intitial condition
 ***********************************/
void solver_static_schrodinger_eq(wavefunction &w, hamiltonian &h, double kx ) 
{
	
	double fn_m;
	double fn;
	double fn_p;		
	
	
	double E0 =  kx*kx/2.;
	
	if (kx>=0){ //Positive Momentum
		for (int i=1; i<w.n1-1; i++) {
			fn_m	    = 2.*h.m1*(E0-h.v[i-1]);
			fn			= 2.*h.m1*(E0-h.v[i]);
			fn_p	    = 2.*h.m1*(E0-h.v[i+1]);
			
			
			complex phim = complex(w.w[i-1][0], w.w[i-1][1] );
			complex phin = complex(w.w[i][0]  , w.w[i][1]   );	
			
			complex phip = (  2.*(1.   -  5.*w.dx1*w.dx1*fn/12. )*phin 
							- ( 1.   +  w.dx1*w.dx1*fn_m/12. )*phim  )
			/( 1.   +  w.dx1*w.dx1*fn_p/12.  );
			
			
			w.w[i+1][0] =	real(phip);
			w.w[i+1][1] =	imag(phip);
			
		}	
	}else { //Negative Momentum
		for (int i=w.n1-2; i>0; i--) {
			fn_m	    = 2.*h.m1*(E0-h.v[i-1]);
			fn			= 2.*h.m1*(E0-h.v[i]);
			fn_p	    = 2.*h.m1*(E0-h.v[i+1]);
			
			
			complex phip = complex(w.w[i+1][0], w.w[i+1][1] );
			complex phin = complex(w.w[i][0]  , w.w[i][1]   );	
			
			complex phim = (  2.*( 1.  -  5.*w.dx1*w.dx1*fn/12. )*phin 
							- ( 1.  +  w.dx1*w.dx1*fn_p/12. )*phip  )
			/( 1.  +  w.dx1*w.dx1*fn_m/12.  );
			
			
			w.w[i-1][0] =	real(phim);
			w.w[i-1][1] =	imag(phim);
			
		}	
	}
	
} //End Numerov TwoWay






//*********************************************************************//
//Building Scattering Continuum wavefunction of momentum k
//******************************************************//
void Continuum_WF( wavefunction &w, hamiltonian &h, double kx )
{
	int n;
	complex psin_1;
	complex Dpsin_1;
	
	//w.expected_kin=kx*kx/2.;	//Kinetic Energy
	
	
	
	if (kx>=0){//POSITIVE MOMENTUM
		n=w.n1-2;
		w.w[0][0] = cos( kx*w.x1[0] );
		w.w[0][1] = sin( kx*w.x1[0] );
		
		w.w[1][0] = cos( kx*w.x1[1] );
		w.w[1][1] = sin( kx*w.x1[1] );			
	}
	else{//NEGATIVE MOMENTUM		
		w.w[w.n1-1][0] = cos( kx*w.x1[w.n1-1] );
		w.w[w.n1-1][1] = sin( kx*w.x1[w.n1-1] );
		
		w.w[w.n1-2][0] = cos( kx*w.x1[w.n1-2] );
		w.w[w.n1-2][1] = sin( kx*w.x1[w.n1-2] );						
		
		
		n=1;		
	}
	
	solver_static_schrodinger_eq(w,h,kx);
	
	psin_1			= complex( w.w[n][0], w.w[n][1] );	
	Dpsin_1			= complex(  diff_th( w.w[n+1][0], w.w[n-1][0], w.dx1 ) 
							  , diff_th( w.w[n+1][1], w.w[n-1][1], w.dx1 )  );
	
	complex factor	= norm_fac( psin_1, Dpsin_1, kx, w.x1[n]);	
	cproduct(w,factor);	
	
	
}




inline double diff_th( double yn1, double yn_1, double h){
	return (yn1-yn_1)/2./h;
}



inline double diff(double yn, double yn_1, double h){
	return (yn-yn_1)/h;
}



inline void cproduct( wavefunction &w, complex a){
	
	for (int i=0; i<w.n1*w.n2*w.n3; i++) {
		complex aux = complex( w.w[i][0], w.w[i][1] );
		aux         = a*aux;
		
		w.w[i][0]	= real(aux);
		w.w[i][1]	= imag(aux);
	}
}



inline complex norm_fac( complex psi, complex Dpsi, double k, double xn)
{
	return 2.*k*complex( cos(k*xn) , sin(k*xn) )/( k*psi - I*Dpsi)/sqrt(dospi);
}






void Coulomb_WF(wavefunction &w, hamiltonian &h, double kx)
{	
	double argument;
	double soft_core=2.;
	for (int i=0; i<w.n1; i++) {
		argument  =  kx*w.x1[i] 
					+ 1./kx*log( abs( w.x1[i]*w.x1[i]
									   + sqrt(w.x1[i]*w.x1[i] +soft_core) ) );
		w.w[i][0] = cos( argument );
		w.w[i][1] = sin( argument );
	}
		
}




//Put dense wdense wave on wundense wavefunction //
void place_dense_WF(wavefunction &wundense, wavefunction &wdense, int skiper1, int skiper2, int skiper3)
{
	
	for(int k=0;k<wundense.n3;k++)
		for(int j=0;j<wundense.n2;j++)
			for(int i=0;i<wundense.n1;i++)
			{
				wundense.w[wundense.index(k,j,i)][0] = wdense.w[wdense.index(k*skiper3,j*skiper2,i*skiper1)][0];
				wundense.w[wundense.index(k,j,i)][1] = wdense.w[wdense.index(k*skiper3,j*skiper2,i*skiper1)][1];
				
			}
}




/** Rotate 90º function ***/
void Rot90(vector<double>& v )
{
	int N=v.size();
	vector<double> a;
	a.resize(N,0.);
	
	for (int i=0; i<N; i++)
		a[i]=v[N-1-i];
	for (int i=0; i<N; i++)
		v[i]=a[i];		
}
/* End rotate 90º */




//Start function mayor
double qmajor(vector<double>& v)
{
	int n=v.size();
	double may = v[0];
	
	for (int k=1;k<n;k++)   
		if (v[k]>may) may=v[k];
	
	return may;
}//End mayor





//Start function minus
double qminus(vector<double>& v)
{
	int n=v.size();
	double min=v[0];
	
	for (int k=1;k<n;k++)
		if (v[k]<min) min=v[k];	
	
	return min;
}//End minus




//Start hmax_ to x
double hmax_x1D(vector<double>& x, vector<double>& v)
{
	int n=v.size();
	double may = v[0];
	double xmax=x[0];
	
	for (int k=1;k<n;k++) 
		if (v[k]>may) {
			may  = v[k];
			xmax = x[k];
		}
	
	return xmax;
}//hmax_ mayor



//Start FWHM to x
double FWHM(vector<double>& x, vector<double>& v)
{
	int N=v.size();
	double hmay  = qmajor(v);
	hmay=hmay/2.0;
	
	double xmax0 = 0. ;
	double xmax1 = 0. ;	
	int flag=0;
	
	for (int i=0; i<N; i++) {
		if (v[i]>=hmay && flag==0) {
			xmax0=x[i];
			flag=1;
		}
		if (v[i]<=hmay && flag==1) {
			xmax1=x[i];
			break;//flag=2;
		}		
	}
	
	return abs(xmax1-xmax0);
}//End FWHM


int kmajor(vector<double>& v)
{
	int n=v.size();
	int kmay=0;  
	double may = qmajor(v);
	
	for (int k=0;k<n;k++)   
		if (v[k]==may)           
			kmay=k;
	
	return kmay;
	
}//End mayor

//Start function minus
int kminus(vector<double>& v)
{
	int n=v.size();
	int kmin=0;
	double min=qminus(v);
	
	for (int k=0;k<n;k++)
		if (v[k]==min)
			kmin=k;
	
	return kmin;
}//End minus




//Tranform unit of atomic unit to SI 
double length_au_SI (double _leng0)
{     
	double leng0=_leng0;     
	double leng1=BohrRadius*leng0;
	
	return leng1;
}

double frequency_length_SI (double _w0)  //from frequency in atomic unit to wavelength International Sistem transformation 
{
	double w0=_w0;
	double lambda=lightC_au*dospi/w0;
	double leng=length_au_SI(lambda);
	
	return leng;
}

double frequency_au_SI (double _w0)  //from frequency in atomic unit to frequency in SI 
{
	double w0=_w0;
	double period_au=dospi/w0;
	
	double period_SI=time_au_SI(period_au);
	double w=dospi/period_SI;
	
	return w;
}

double time_au_SI (double _time0)
{
	double time0=_time0;
	double time1=time_SI*time0;
	
	return time1;     
}

double intensity_au_SI (double _I0)
{
	double I0=_I0;
	double I1=I0*3.5e16;
	
	return I1;     
}

//Tranform unit of SI to atomic unit 
double length_SI_au (double _leng0)
{     
	double leng0=_leng0;     
	double leng1=leng0/BohrRadius;
	
	return leng1;
}

//Wavelength of SI unit to atomic unit frequency
double length_frequency_au (double _lambda)  
{
	double lambda_SI = _lambda;
	double lambda_au = length_SI_au(lambda_SI);
	double w0 = lightC_au*dospi / lambda_au;
	
	return w0;
}

double frequency_SI_au (double _w0)  
{
	double w0 = _w0;
	double period_SI = dospi/w0;
	double w = dospi / time_SI_au (period_SI);
	
	return w;
}

double time_SI_au (double _time0)
{
	double time0=_time0;
	double time1=time0/time_SI;
	
	return time1;     
}

double intensity_SI_au (double _I0)
{
	double I0=_I0;
	double I1=I0/3.5e16;
	
	return I1;     
}

