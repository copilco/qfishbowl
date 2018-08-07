#ifndef TIMEOBJECT_H
#define TIMEOBJECT_H
#include <iostream>
#include "fftw3.h"
#include<vector>
//!  A Hamiltonian class. 
/*!
 Here all the variables and methods for the 
wavefunction are specified
 */


class timeobject: public timegrid
{
public:
	
	int space_indicator;
	
	fftw_complex *f;
	fftw_plan pF, pB;

	
	double wnorm;
	double expected_kin;
	
	int something;
	double fR_double;
	double fI_double;
	
	void initialize()
    {
		f = (fftw_complex*) fftw_malloc(n * sizeof(fftw_complex)); //reserva espacio de memoria
		if(!f) cout<<"\nError allocating array, please check ...";
		
		pF  =fftw_plan_dft_1d(n,f,f, 1, FFTW_ESTIMATE);
		pB  =fftw_plan_dft_1d(n,f,f,-1, FFTW_ESTIMATE);
		
		for(int i=0;i<n;i++)
		{
			f[i][0]=0.; 
			f[i][1]=0.;
		}
		
		wnorm=0.;
		// expected_kin=0.;
		fR_double=0.0;
		fI_double=0.0;
		
		//p1DF_1 = fftw_plan_dft_1d(n1,f[0], f[0], FFTW_FORWARD, FFTW_ESTIMATE);
    }
	
	
	void put_on_grid(timegrid &g)
    {
		set_grid(g.n,g.dt,g.t0 );
		initialize();
    }

	
	/******************************************************/
	// Functions of the wavefunction itself
	double norm()
    {
      /* Norm */
		double norm0=0.0;
		
		for(int i=0;i<n;i++)
			norm0+=dt*(f[i][0]*f[i][0]+f[i][1]*f[i][1]);
		return norm0;  
    }
	
	/******************************************************/
	//Normalize.
	void normalize()
    {
		double norm0=norm(); //que se entiende en este paso?
		//cout << "\n norm imput "<<norm0;
		for(int i=0;i<n;i++)
	{
		f[i][0]/=sqrt(norm0);
		f[i][1]/=sqrt(norm0);
	}
		
		norm0=norm();
		//cout << "\n norm output "<<norm0;
		
    }
	
	/******************************************************/
	//Expected value x1
	double expected_t()
    {
		double expec_t=0.;
		
		for(int i=0;i<n;i++)
	    { 
	      expec_t+=dt*(f[i][0]*t[i]*f[i][0]+f[i][1]*t[i]*f[i][1]);
	    }
		
		return expec_t;
    }
	
	
	double w_norm()
	{

		wnorm=0.;
		double factor=1./ n;
		
		/**********************/
		fftw_execute( pF);
		/**********************/
		
		for(int i=0;i< n;i++)
		{
			wnorm+= dw*
			( f[ i][0]* f[ i][0]+ f[ i][1]* f[ i][1])
			* kfactor* kfactor;
			
			f[i][0]*=factor;
			f[i][1]*=factor;
	    
		}
		
		/**********************/
		fftw_execute( pB);
		/**********************/
		return wnorm;
  }
	
	void integrate()
	{
		
		fftw_execute(pF);
		//		double dw=2.*pi/nmaxt/dt;
		
		
		f[0][0]=0.;
		f[0][1]=0.;
		
		for(int ktime=1;ktime<n;ktime++)
		{
			
			double tem0=f[ktime][0];
			double tem1=f[ktime][1];
			
			f[ktime][0]= tem1/w[ktime]/n;
			f[ktime][1]=-tem0/w[ktime]/n;
		}
		fftw_execute(pB);
		double dc=(f[0][0]+f[0][1])/n;
		
		for(int ktime=1;ktime<n;ktime++)
		{
			
			double tem0=f[ktime][0];
			double tem1=f[ktime][1];
			
			f[ktime][0]= tem0+dc*t[ktime];
			f[ktime][1]= tem1+dc*t[ktime];
		}
		
		
	}

	
	/*** RungeKuttaIntegrator  ***/

	void integrateRK4(){
		//Double time grid and electric field
		complex alpha0 = complex(0.0,0.0);
		complex alpha1 = complex(0.0,0.0);
		complex alpha2 = complex(0.0,0.0);
		complex alpha3 = complex(0.0,0.0);		
		
		
		vector<complex> v(n,0.0);
		
		for (int ktime=0; ktime<n; ktime++)
			v[ktime]=complex(f[ktime][0],f[ktime][1]);
		
		
		f[0][0]=0.0;
		f[0][1]=0.0;
		
		//Implementation RK4th		
		for (int ktime=0; ktime<n-1 ; ktime++) {
			
			alpha0 = v[ktime];
			 
			
			fR_double=real(v[ktime]) + ( real(v[ktime+1]) -real(v[ktime]) )/2.0 ;	
			fI_double=imag(v[ktime]) + ( imag(v[ktime+1]) -imag(v[ktime]) )/2.0 ;	
			
			
			alpha1 = complex(fR_double,fI_double);
			
			
			alpha2 = alpha1;
			
			
			alpha3 = v[ktime+1];	
			
			
			f[ktime+1][0] = ( f[ktime][0] + dt/6.0*(real(alpha0 + 2.0*(alpha1 + alpha2) + alpha3) ) );

			
			f[ktime+1][1] = ( f[ktime][1] + dt/6.0*(imag(alpha0 + 2.0*(alpha1 + alpha2) + alpha3) ) );		
		
		}
			
		
	}

	
	void differentiate()
	{
		fftw_execute(pF);
		
		for(int ktime=1;ktime<n;ktime++)
		{
			
			double tem0=f[ktime][0];
			double tem1=f[ktime][1];
			
			f[ktime][0]=-tem1*w[ktime]/n;
			f[ktime][1]= tem0/w[ktime]/n;
		}
		
		fftw_execute(pB);
		
	}
	
};
#endif
