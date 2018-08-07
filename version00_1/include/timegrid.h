// *** grid.h ***  for the hop_cart project -----------------  db 09/98
#ifndef TIMEGRID_H
#define TIMEGRID_H
#include<iostream>
#include<math.h>
#include<vector> 
#include <fftw3.h>
#include "constant.h"
using namespace std;

//!  A timegrid class. 
/*!
 This class wil be used for the creation of pulses and quantities related to
 evolution in time. It will contain
 */

class timegrid
	{
	public:
		
		//string name;
		int n;   //*!< Number of points in the x1 direction */ 
		
		double dt;  /*!< Value of the dx in x1 */ 
		double t0;
		
		double nisq; /*!< Maximum momentum resolved in the x1 direction Nysquit freq.*/
		double dw;   /*!< Delta of momentum in the q1 direction.*/		
		
		double kfactor; /*!< Factor for the right value after the fft on x1. It is tricky but it works well*/            
		
		vector<double> t; /*!< The vector containing the x1 coordinate frame */
		vector<double> w; /*!< The vector containing the w coordinate frame */		
		
		/********************************************/
		// Useful data
		
		inline double sizet_au() const { return dt*n; }
		inline double sizew_au() const { return dw*n; }
		
		inline double vol_elem() const { return dt; }
		inline double qvol_elem() const { return dw; }
		
		inline long tot_pts() const { return n; }
		inline double aprox_size_Mb() const{ return n*16.*1e-6; }
		
		
		/********************************************/
		// Initialize
		inline void set_dt(double _t0, double _dt)
		  {
			dt=_dt;
			t0=_t0;
		  }
		
		inline void set_n(int _n)
		  {
			n=_n;
		  }
		
		inline void set_grid(int _n, double _dt, double _t0)
		  {			  
			  set_n(  _n);
			  set_dt(_t0, _dt);
			
			  t.resize(n, 0.);			
			  w.resize(n, 0.);
			  
			  nisq=pi/dt;
			  dw=dospi/n/dt;
			  kfactor=dt/sqrt(dospi);

			  for(int i=0;i<n;i++)
				  t[i] = t0 + i*dt;	//t[i]	=-n*dt/2+i*dt;
			  
			  /**************************/
			  for(int i=0;i<=n/2;i++)
				  w[i]=i*dw;
			
			  for(int i=n/2+1;i<n;i++)
				  w[i]=(i-n)*dw;//-nisq+(i-n/2)*dw;
		  }	
}; //End of the grid.h 
#endif
