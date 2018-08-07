#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <iostream>
#include "fftw3.h"


class wavefunction: public grid
{
	
public:
	
	int space_indicator;

	//int nthreads;	
	//int fftw_init_threads(void);

	
	
	fftw_complex *w;
	fftw_plan p3DF, p3DB;
	
	double qnorm;
	double expected_kin;
	double pmom1;	
	double pmom2;
	double pmom3;
	double norm_factor;
	int something;
	
	vector<double> phase;
	vector<double> Qphase;
	vector<double> falter;
	//  vector<double> prob;
	
	
	
	void cleanObjects(){
		if (w!=NULL)	fftw_free(w);			
	}
	
	
	
	/*~wavefunction(){
		cleanObjects();
	}
	//*/
	
	
	void initialize()
	{
		
		//fftw_plan_with_nthreads(nthreads);	
		w = (fftw_complex*) fftw_malloc(n1*n2*n3 * sizeof(fftw_complex)); ///malloc it is a change of size of w
		if(!w) cout<<"\nError allocating array, please check ...";
		
		
		phase.resize(n1*n2*n3,	0.);
		Qphase.resize(n1*n2*n3,	0.);
		falter.resize(n1*n2*n3,	1.);
		
		p3DF  =fftw_plan_dft_3d(n3,n2,n1,w,w,-1, FFTW_ESTIMATE);  //Directy Fourier Transform 
		p3DB  =fftw_plan_dft_3d(n3,n2,n1,w,w, 1, FFTW_ESTIMATE);  //Inverse Fourier Transform    

		
		qnorm			=0.;
		expected_kin	=0.;
		alter();
		norm_factor		=1.;
		
	}
	
	
	void put_on_grid(grid &g)
	{
		set_grid(g.n1,g.n2,g.n3,g.dx1,g.dx2,g.dx3 );
		
		symmetry_n1=g.symmetry_n1;
		symmetry_n2=g.symmetry_n2;
		symmetry_n3=g.symmetry_n3;
		
		initialize();
	}
	
	
	//Direct Fourier Transform 
	void Direct_Fourier_Transform()
	{
		
		/**********************/
		fftw_execute( p3DF);
		/**********************/    

	}
	
	
	//Back Fourier Transform 
	void Back_Fourier_Transform()
	{		

		/**********************/
		fftw_execute( p3DB);
		/**********************/
	}	
	
	
	
	void alter()
	{
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{
					int a=index(k,j,i);
					if(a%2!=0)
						falter[a]*=-1;		
				}
	}
	
	

	
	/******************************************************/
	// Functions of the wavefunction itself
	double norm() 
	{
		double norm0=0.0;      
		for(int i=0;i<n1*n2*n3;i++)
			norm0+=dx1*dx2*dx3*(w[i][0]*w[i][0]+w[i][1]*w[i][1]);      
		return norm0;
	}
	
	//q_norm
	double q_norm()
	{
		qnorm=0.;
		double factor=1./ n1/ n2/ n3;    
		/**********************/
		fftw_execute( p3DF);
		/**********************/    
		for(int k=0;k< n3;k++)
			for(int j=0;j< n2;j++)
				for(int i=0;i< n1;i++)
				{
					qnorm+= dq1* dq2* dq3*
					( w[ index(k,j,i)][0]* w[ index(k,j,i)][0]+ w[ index(k,j,i)][1]* w[ index(k,j,i)][1])
					* kfactor1* kfactor1* kfactor2* kfactor2* kfactor3* kfactor3;
					
					w[ index(k,j,i)][0]*=factor;
					w[ index(k,j,i)][1]*=factor;	    
				}    
		/**********************/
		fftw_execute( p3DB);
		/**********************/
		return qnorm;
	}
	
	
	
	double qq_norm()
	{
		qnorm=0.;    
		for(int i=0;i< n1*n2*n3;i++)
		{
			qnorm+= dq1* dq2* dq3*
			( w[i][0]* w[i][0] + w[i][1]* w[i][1]);
		}    
		return qnorm;
	}
	
	
	
	
	
	double wq_norm()
	{
		qnorm=0.;    
		double factor=1./ n1/ n2/ n3;    
		
		fftw_execute( p3DF);
		for(int i=0;i< n1*n2*n3;i++)
		{
			qnorm+= dq1* dq2* dq3*
			( w[i][0]* w[i][0] + w[i][1]* w[i][1]);
			w[i][0]*=factor;
			w[i][1]*=factor;
		}
		fftw_execute( p3DB);
		
		return qnorm;
	}

	
	
	//Normalize.
	void normalize()
	{
		double norm0=norm()/norm_factor;
		for(int i=0;i<n1*n2*n3;i++)
		{
			w[i][0]/=sqrt(abs(norm0));
			w[i][1]/=sqrt(abs(norm0));
		}
	}

	
	//STRAING FUNCTION
	void normalize_q()
	{
		double qnorm0=q_norm();
		double qnorm1=wq_norm();
		double factor=1./ n1/ n2/ n3;
		
		/**********************/
		fftw_execute( p3DF);
		/**********************/
		
		for(int i=0;i<n1*n2*n3;i++)
		{
			w[i][0]/=sqrt(qnorm0);
			w[i][1]/=sqrt(qnorm0);
			

			
			w[i][0]*=factor;
			w[i][1]*=factor;
		}      
		
		/**********************/
		fftw_execute( p3DB);
		/**********************/
	}

	
	
	void normalize_qq()
	{
		double qnorm0=qq_norm();
		for(int i=0;i<n1*n2*n3;i++)
		{
			w[i][0]/=sqrt(qnorm0);
			w[i][1]/=sqrt(qnorm0);
		}         
	}

	
	
	
	
	
	/******************************************************/
	//Expected value x1
	double expected_x1()
    {
		double expec_x1=0.;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					expec_x1+=dx1*dx2*dx3*(w[index(k,j,i)][0]*x1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*x1[i]*w[index(k,j,i)][1]);
				}
		
		return expec_x1;
    }

	
	
	
	double expected_halfx1(int sw0)
    {
		double expec_x1=0.;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					if (i>=n1/2 && sw0==1)
						expec_x1+=dx1*dx2*dx3*(w[index(k,j,i)][0]*x1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*x1[i]*w[index(k,j,i)][1]);
					
					if (i<n1/2 && sw0==0)
						expec_x1+=dx1*dx2*dx3*(w[index(k,j,i)][0]*x1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*x1[i]*w[index(k,j,i)][1]);
				}
		
		return expec_x1;
    }  

	
	
	
	double expected_absx1()
    {
		double expec_x1=0.;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					expec_x1+=dx1*dx2*dx3*(w[index(k,j,i)][0]*abs(x1[i])*w[index(k,j,i)][0]+w[index(k,j,i)][1]*abs(x1[i])*w[index(k,j,i)][1]);
				}
		
		return expec_x1;
    }
	
	
	
	
	//Expected value x2
	double expected_x2()
    {
		double expec_x2=0.;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					expec_x2+=dx1*dx2*dx3*(w[index(k,j,i)][0]*x2[j]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*x2[j]*w[index(k,j,i)][1]);
				}
		
		return expec_x2;
    }

	
	
	
	//Expected value x3
	double expected_x3()
    {
		double expec_x3=0.;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					expec_x3+=dx1*dx2*dx3*(w[index(k,j,i)][0]*x3[k]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*x3[k]*w[index(k,j,i)][1]);
				}
		
		return expec_x3;
    }

	
	
	
	
	/******************************************************/
	//Expected value q1
	double expected_q1() 
    {
		/*************************/
		fftw_execute( p3DF);
		/************************/  
		double expec_q1=0.;
		double qq_norm = q_norm();
		double factor=1./ n1/ n2/ n3;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					expec_q1+=dq1*dq2*dq3*(w[index(k,j,i)][0]*q1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q1[i]*w[index(k,j,i)][1])		
					*kfactor1*kfactor1*kfactor2*kfactor2*kfactor3*kfactor3;
					
					w[ index(k,j,i)][0]*=factor;
					w[ index(k,j,i)][1]*=factor;	    
				}
		
		/**********************/
		fftw_execute( p3DB);
		/**********************/
		return expec_q1;
		
    }//End expected value q1

	
	
	
	//Expectedq value q1
	double expectedq_q1() 
    {
		double expec_q1=0.;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					expec_q1+=dq1*dq2*dq3*(w[index(k,j,i)][0]*q1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q1[i]*w[index(k,j,i)][1]);
				}
		
		return expec_q1;
		
    }//End expectedq value q1
	
	
	
	
	
	//Expected value q1 half
	double expected_halfq1(int sw0) 
    {
		/*************************/
		fftw_execute( p3DF);
		/************************/  
		double expec_q1=0.;
		double qq_norm = q_norm();
		double factor=1./ n1/ n2/ n3;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					if ( i<n1/2 && sw0==1 )		
						expec_q1+=2*dq1*dq2*dq3*(w[index(k,j,i)][0]*q1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q1[i]*w[index(k,j,i)][1])		
						*kfactor1*kfactor1*kfactor2*kfactor2*kfactor3*kfactor3;
					
					if ( i>=n1/2 && sw0==0 )		
						expec_q1+=2*dq1*dq2*dq3*(w[index(k,j,i)][0]*q1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q1[i]*w[index(k,j,i)][1])		
						*kfactor1*kfactor1*kfactor2*kfactor2*kfactor3*kfactor3;
					
					w[ index(k,j,i)][0]*=factor;
					w[ index(k,j,i)][1]*=factor;	    
				}
		
		/**********************/
		fftw_execute( p3DB);
		/**********************/
		return expec_q1;
		
    }//End expected value q1 half

	
	
	
	
	//Expectedq value q1 half
	double expectedq_halfq1(int sw0) 
    {
		double expec_q1=0.;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					if (i<n1/2 && sw0==1)		
						expec_q1+=2*dq1*dq2*dq3*(w[index(k,j,i)][0]*q1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q1[i]*w[index(k,j,i)][1]);
					
					if (i>=n1/2 && sw0==0)		
						expec_q1+=2*dq1*dq2*dq3*(w[index(k,j,i)][0]*q1[i]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q1[i]*w[index(k,j,i)][1]);
				}
		
		return expec_q1;
    }//End expectedq value q1 half

	
	

	
	
	
	//Expected value q2
	double expected_q2()
    {
		/*************************/
		fftw_execute( p3DF);
		/************************/  
		double expec_q2=0.;
		double factor=1./ n1/ n2/ n3;
		
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					expec_q2+=dq1*dq2*dq3*(w[index(k,j,i)][0]*q2[j]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q2[j]*w[index(k,j,i)][1])		
					*kfactor1*kfactor1*kfactor2*kfactor2*kfactor3*kfactor3;
					
					w[ index(k,j,i)][0]*=factor;
					w[ index(k,j,i)][1]*=factor;	    
				}
		
		/**********************/
		fftw_execute( p3DB);
		/**********************/
		return expec_q2;
		
    }//End expected value q2

	

	
	
	
	//Expected value q2 half
	double expected_halfq2(int sw0) 
    {
		/*************************/
		fftw_execute( p3DF);
		/************************/  
		double expec_q2=0.;
		double qq_norm = q_norm();
		double factor=1./ n1/ n2/ n3;
		
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					if (j<n2/2 && sw0==1)
					{
						expec_q2+=2*dq1*dq2*dq3*(w[index(k,j,i)][0]*q2[j]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q2[j]*w[index(k,j,i)][1])		
						*kfactor1*kfactor1*kfactor2*kfactor2*kfactor3*kfactor3;
					}
					if (j>=n2/2 && sw0==0)
					{
						expec_q2+=2*dq1*dq2*dq3*(w[index(k,j,i)][0]*q2[j]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q2[j]*w[index(k,j,i)][1])		
						*kfactor1*kfactor1*kfactor2*kfactor2*kfactor3*kfactor3;
					}
					
					w[ index(k,j,i)][0]*=factor;
					w[ index(k,j,i)][1]*=factor;	    
				}
		
		/**********************/
		fftw_execute( p3DB);
		/**********************/
		return expec_q2;
		
    }//End expected value q2

	
	
	//Expected value q3
	double expected_q3()
    {
		/*************************/
		fftw_execute( p3DF);
		/************************/  
		double expec_q3=0.;
		double qq_norm = q_norm();
		double factor=1./ n1/ n2/ n3;
		
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					expec_q3+=dq1*dq2*dq3*(w[index(k,j,i)][0]*q3[k]*w[index(k,j,i)][0] + w[index(k,j,i)][1]*q3[k]*w[index(k,j,i)][1])		
					*kfactor1*kfactor1*kfactor2*kfactor2*kfactor3*kfactor3;
					
					w[ index(k,j,i)][0]*=factor;
					w[ index(k,j,i)][1]*=factor;	    
				}
		
		/**********************/
		fftw_execute( p3DB);
		/**********************/
		return expec_q3;
		
    }//End expected value q3

	
	
	//Expected value q3 half
	double expected_halfq3() 
    {
		/*************************/
		fftw_execute( p3DF);
		/************************/  
		double expec_q3=0.;
		double qq_norm = q_norm();
		double factor=1./ n1/ n2/ n3;
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					if (k<n3/2)
					{
						expec_q3+=2*dq1*dq2*dq3*(w[index(k,j,i)][0]*q3[k]*w[index(k,j,i)][0]+w[index(k,j,i)][1]*q3[k]*w[index(k,j,i)][1])		
						*kfactor1*kfactor1*kfactor2*kfactor2*kfactor3*kfactor3;
					}
					
					w[ index(k,j,i)][0]*=factor;
					w[ index(k,j,i)][1]*=factor;	    
				}
		
		/**********************/
		fftw_execute( p3DB);
		/**********************/
		return expec_q3;
		
    }//End expected value q3 half  

		
	
	
	//=====================================//
	//Definition of the zero on the wavefunction 
	//position and momentum space
	void wzero(double zeroR, double zeroI)
	{
		for(int i=0;i<n1*n2*n3;i++)
		{
            if(abs(w[i][0])<=zeroR) w[i][0]=0.0;		
            if(abs(w[i][1])<=zeroI) w[i][1]=0.0;		
		}
	}
	

	
	double atan2Mod(double y, double x, double zeroY, double zeroX)
	{
		if(abs(x)<=zeroX) x=0.0;
		if(abs(y)<=zeroY) y=0.0;
		return atan2(y,x);
	}
	
	

};

#endif
