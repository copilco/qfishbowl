#include <iostream>
//#include <grid.h>
#include "fftw3.h"



class wavefunction: public grid
{

 public:
  
  int space_indicator;

  fftw_complex *w;
  fftw_plan p3DF, p3DB;
  
  double qnorm;
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
      w = (fftw_complex*) fftw_malloc(n1*n2*n3 * sizeof(fftw_complex));
      if(!w) cout<<"\nError allocating array, please check ...";
      
      p3DF  =fftw_plan_dft_3d(n3,n2,n1,w,w,-1, FFTW_ESTIMATE);
      p3DB  =fftw_plan_dft_3d(n3,n2,n1,w,w, 1, FFTW_ESTIMATE);

      qnorm=0.;
      expected_kin=0.;

      //p1DF_1 = fftw_plan_dft_1d(n1,w[0], w[0], FFTW_FORWARD, FFTW_ESTIMATE);
    }

  void put_on_grid(grid g)
    {

      set_grid(g.n1,g.n2,g.n3,g.dx1,g.dx2,g.dx3 );
          
      symmetry_n1=g.symmetry_n1;
      symmetry_n2=g.symmetry_n2;
      symmetry_n3=g.symmetry_n3;
      
      initialize();
    }

  /******************************************************/
  // Functions of the wavefunction itself
  double norm()
    {
      /* Norm */
      double norm0=0.0;
      
      for(int i=0;i<n1*n2*n3;i++)
	norm0+=dx1*dx2*dx3*(w[i][0]*w[i][0]+w[i][1]*w[i][1]);
      
      return norm0;  
    }

  /******************************************************/
  //Normalize.

  void normalize()
    {
      double norm0=norm();
      //cout << "\n norm imput "<<norm0;
      for(int i=0;i<n1*n2*n3;i++)
	{
	  w[i][0]/=sqrt(norm0);
	  w[i][1]/=sqrt(norm0);
	}
      
      norm0=norm();
      //cout << "\n norm output "<<norm0;
     
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
  
  //Expected value x1
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
  
  

};
