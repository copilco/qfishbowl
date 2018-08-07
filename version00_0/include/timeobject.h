#include <iostream>
//#include <grid.h>
#include "fftw3.h"


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


  /*
  fftw_plan p1DF_1;      
  fftw_plan p1DB_1;      

  fftw_plan p1DF_2;      
  fftw_plan p1DB_2;      
  
  fftw_plan p1DF_3;      
  fftw_plan p1DB_3;      
  */
  int something;
  
  
  void initialize()
    {
		f = (fftw_complex*) fftw_malloc(n * sizeof(fftw_complex));
		if(!f) cout<<"\nError allocating array, please check ...";
		
		pF  =fftw_plan_dft_1d(n,f,f,-1, FFTW_ESTIMATE);
		pB  =fftw_plan_dft_1d(n,f,f, 1, FFTW_ESTIMATE);
		
		for(int i=0;i<n;i++)
		{
			f[i][0]=0.;
			f[i][1]=0.;
		}
		
      wnorm=0.;
     // expected_kin=0.;

      //p1DF_1 = fftw_plan_dft_1d(n1,f[0], f[0], FFTW_FORWARD, FFTW_ESTIMATE);
    }

  void put_on_grid(timegrid g)
    {
		set_grid(g.n,g.dt );
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
      double norm0=norm();
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
