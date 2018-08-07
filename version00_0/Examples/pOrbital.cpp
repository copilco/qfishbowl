//#include <iostream>
#include <iostream>
#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "cartesian.h"
//#include "constant.h"
using namespace std;
//#include <complex.h>
//#define complex complex<double>

int main()
{
  
  /****************************/
  //These are grid parameters

  int nx=30;
  int ny=30;
  int nz=30;
  
  double dx=0.3;
  double dy=0.3;
  double dz=0.3;//1.;
  /****************************/

  grid g1;                                    //Define a grid object
  g1.set_grid(nx,ny,nz,dx,dy,dz);             //Inittialize  
  
  //Sone info

  cout << "Size in n1  "<< g1.n1 << "\n";     
  cout << "Vol elemnt  "<< g1.vol_elem() << "\n";
  cout << "Size Mb     "<< g1.aprox_size_Mb() << "\n";

  /****************************/  
  //A wavefunction can be initialized using a grid as parameter.
  
  wavefunction w1,w2; 
  w1.put_on_grid(g1);
  w2.put_on_grid(g1);

  hamiltonian h;
  h.put_on_grid(g1);
  h.set_v_hlike(1., 0.6 );

  /*****************************************************/ 
  // Start the gaussian
  gaussian(w1);
  gaussian(w2);
  
  w1.normalize(); //Normalize the wavefunction
  w2.normalize(); //Normalize the wavefunction
  
  /*****************************************************/
  //Do the FFT in 3D
  
  FILE *out0,*out1,*out2,*out3,*out4,*out5,*out6,*out7;

  out0=fopen("out0.txt","w");
  out1=fopen("out1.txt","w");
  out2=fopen("out2.txt","w");
  out3=fopen("out3.txt","w");  
  out4=fopen("out4.txt","w");
  out5=fopen("out5.txt","w");
  out6=fopen("out6.txt","w");
  out7=fopen("out7.txt","w");

 
 /*************************************************************************************************/
  complex dt=0.05;//-I*0.05;
  
  
  double e0=sqrt(2e14/3.5e16);
  double w=0.057;
  double period=dospi/w;
  double cycle_number=6.;
  double pol_angle=6.*pi/6.;

  int nmaxt=int(period*cycle_number/abs(dt));
  int snaper1=int(0.125*period/abs(dt));

  cout<<"\n nmaxt "<<nmaxt;
  
  double avect_x=0.;
  double avect_y=0.;

  int snap=0;
  cout<<"\n nmaxt "<<nmaxt;
  /*************************************************************************************************/
  int skiper1=5;
  int skiper2=5;
  int skiper3=1;

  double obsk;
  
  for(int hh=0;hh<=0;hh++)
    {
      if(hh==0)
	{
	  dt=-I*0.05;
	  nmaxt=1000;
	}
      if(hh==1)
	{
	  dt=0.05;
	  nmaxt=int(period*cycle_number/abs(dt));
	}
      if(hh==2)
	dt=0.05;
      
      for(int ktime=0;ktime<=nmaxt;ktime++)
	{

	  
	  complex t = ktime*1.*dt;
	  complex tt=(ktime+0.5)*1.*dt;
	  
	  if(hh==1)
	    {
	      double efield_x=e0*cos(pol_angle)*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t));
	      double efield_y=e0*sin(pol_angle)*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t));
	      
	      avect_x+=-efield_x*abs(dt)*lightC_au;
	      avect_y+=-efield_y*abs(dt)*lightC_au;
	      
	      //	      fprintf(fpt3,"%e %e %e\n",abs(t),efield_x,avect_x/lightC_au);
	      fprintf(out3,"%e %e %e %e %e\n",abs(t),efield_x,avect_x/lightC_au,efield_y,avect_y/lightC_au);
	    }
	  
	  
	  if(hh==0)//Imaginary time
	    {
	      obsk=prop_kinetic(w1, h, dt/2.);
	      //    prop_potential(w1,h,dt);
	      obsk=prop_kinetic(w1, h, dt/2.);
	      // w1.normalize();

	      obsk=prop_kinetic(w2, h, dt/2.);
	      prop_potential(w2,h,dt);
	      obsk=prop_kinetic(w2, h, dt/2.);
	      w2.normalize();

	      for(int k=0;k<w2.n3;k++)
	      for(int j=0;j<w2.n2/2;j++)
		for(int i=0;i<w2.n1;i++)
		  {
		    w2.w[w2.index(k,j,i)][0]=-(w2.w[w2.index(k,(w2.n2-1)-j,i)][0]);
		    w2.w[w2.index(k,j,i)][1]=-(w2.w[w2.index(k,(w2.n2-1)-j,i)][1]);
		  }

	    }//end hh==0
      /*******************************************************
       if((ktime%snaper1)==0)
	 {
	   snap+=1;
	   cout << "snap number " << snap;
	   for(int k=0;k<w1.n3/skiper3;k++)
	     for(int j=0;j<w1.n2/skiper2;j++)
	       for(int i=0;i<w1.n1/skiper1;i++)
		 {
		   
		   double norm=
		     w1.dx1*w1.dx2*w1.dx3*(w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]+
					   w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]);
		   fprintf(out5,"%e\n",norm);
		 }
	 }
       /*******************************************************/
       
       fprintf(out6,"%e %e %e\n",abs(t),w1.expected_x1(),w1.expected_x2());
       cout << "\n Norm after the way on and way back "<<w1.norm();
       cout << "  Kinetic Energy "<< kinetic_finite_diff(w1,h);
       cout << "  Potential energy "<< potential_energy(w1,h);
       cout << "  Energy "<< kinetic_finite_diff(w1,h)+potential_energy(w1,h);
       cout << "  obsk "<< obsk;
      

	}//end time loop
    }

}//End
