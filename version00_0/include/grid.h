// *** grid.h ***  for the hop_cart project -----------------  db 09/98
//#ifndef GRID_H
//#define GRID_H
#include "constant.h"
#include<iostream.h>
#include<vector.h> 
#include "fftw3.h"
using namespace std;

class grid
{

 public:
  
  //string name;
  int n1;
  int n2;
  int n3;

  double dx1;
  double dx2;
  double dx3;
  
  int symmetry_n1;
  int symmetry_n2;
  int symmetry_n3;
   
  double nisq1;
  double dq1;
  
  double nisq2;
  double dq2;

  double nisq3;
  double dq3;

  double kfactor1;
  double kfactor2;
  double kfactor3;

  vector<double> x1;
  vector<double> x2;
  vector<double> x3;

  vector<double> q1;
  vector<double> q2;
  vector<double> q3;
  
  /********************************************/
  //Index in 1D
  
  int index1_1d(int a) const
    {
      int index=0;
      if (a>=0 && a <n1 )
	index=a;
      else
	cout << "Bad index in the 1D n1\n";
      return index;      
    }

  int index2_1d(int a) const
    {
      int index=0;
      if (a>=0 && a <n2 )
	index=a;
      else
	cout << "Bad index in the 1D n2\n";
      return index;      
    }

  int index3_1d(int a) const
    {
      int index=0;
      if (a>=0 && a <n3 )
	index=a;
      else
	cout << "Bad index in the 1D n3\n";
      return index;      
    }
  
  /********************************************/
  //Two dimensional index.

  int index12_2d(int a, int b)
    {
     
      int index=a*n1+b;
      return index;
    }

  int index13_2d(int a, int b)
    {
      int index=a*n2*n1+b;
      return index;
    }

  int index23_2d(int a, int b)
    {
      int index=a*n2*n1+b*n1;
      return index;
    }

  /********************************************/
  // Index for the 3D array

  int index(int a,int b,int c)
    {
      int index=a*n2*n1+b*n1+c;
      return index;
    }
  
  
  /********************************************/
  // Useful data

  inline double size1_au() const { return dx1*n1; }
  inline double size2_au() const { return dx2*n2; }
  inline double size3_au() const { return dx3*n3; }

  inline double sizeq1_au() const { return dq1*n1; }
  inline double sizeq2_au() const { return dq2*n2; }
  inline double sizeq3_au() const { return dq3*n3; }
  
  inline double vol_elem() const { return dx1*dx2*dx3; }
  inline double qvol_elem() const { return dx1*dx2*dx3; }

  inline long tot_pts() const { return n1*n2*n3; }
  inline double aprox_size_Mb() const{ return n1*n2*n3*16.*1e-6; }

  /********************************************/
  // Initialize

  inline void set_dx(double _dx1, double _dx2, double _dx3)
    {
      dx1=_dx1;
      dx2=_dx2;
      dx3=_dx3;
    }

  inline void set_n(int _n1,int _n2,int _n3)
    {
      n1=_n1;
      n2=_n2;
      n3=_n3;
    }

  inline void set_grid(int _n1, int _n2, int _n3, double _dx1, double _dx2, double _dx3)
    {

      set_n(  _n1, _n2, _n3);
      set_dx(_dx1, _dx2,_dx3);
      
      /*
      dx1=_dx1;
      dx2=_dx2;
      dx3=_dx3;

      n1=_n1;
      n2=_n2;
      n3=_n3;      
      */

      x1.resize(n1, 0.);
      x2.resize(n2, 0.);
      x3.resize(n3, 0.);

      q1.resize(n1, 0.);
      q2.resize(n2, 0.);
      q3.resize(n3, 0.);


      if(n1==1)
	{
	  nisq1=pi/dx1;
	  dq1=1.;
	  kfactor1=1.;
	}
      else
	{
	  nisq1=pi/dx1;
	  dq1=dospi/n1/dx1;
	  kfactor1=dx1/sqrt(dospi);
	}
      if(n2==1)
	{
	  nisq2=pi/dx2;
	  dq2=1.;
	  kfactor2=1.;
	}
      else
	{
	  nisq2=pi/dx2;
	  dq2=dospi/n2/dx2;
	  kfactor2=dx2/sqrt(dospi);
	}

      if(n3==1)
	{
	  nisq3=pi/dx3;
	  dq3=1.;
	  kfactor3=1.;
	}
      else
	{
	  nisq3=pi/dx3;
	  dq3=dospi/n3/dx3;
	  kfactor3=dx3/sqrt(dospi);
	}

      /*
      nisq1=pi/dx1;
      dq1=dospi/n1/dx1;
      
      nisq2=pi/dx2;
      dq2=dospi/n2/dx2;
      
      nisq3=pi/dx3;
      dq3=dospi/n3/dx3;
      */

      for(int k=0;k<n3;k++)
	x3[k]=(-(n3)/2.+k)*dx3;
      
      for(int j=0;j<n2;j++)
	x2[j]=(-(n2)/2.+j)*dx2;
      
      for(int i=0;i<n1;i++)
	x1[i]=(-(n1)/2.+i)*dx1;
      
      /**************************/
      
      for(int k=0;k<n3/2;k++)
	q3[k]=(k)*dq3;
      
      for(int k=n3/2;k<n3;k++)
	q3[k]=-nisq3+(k-n3/2)*dq3;
      
      for(int j=0;j<n2/2;j++)
	q2[j]=(j)*dq2;
      
      for(int j=n2/2;j<n2;j++)
	q2[j]=-nisq2+(j-n2/2)*dq2;
      
      for(int i=0;i<n1/2;i++)
	q1[i]=(i)*dq1;
      
      for(int i=n1/2;i<n1;i++)
	q1[i]=-nisq1+(i-n1/2)*dq1;

    }

  //End of the grid.h 
};

