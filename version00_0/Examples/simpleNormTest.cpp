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

  int nx=50;
  int ny=50;
  int nz=50;
  
  double dx=0.15;
  double dy=0.15;
  double dz=0.15;
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
  complex dt=0.05;//-I*0.05;//-I*0.05;
  
  
  double e0=sqrt(2e14/3.5e16);
  double w=0.057;
  double period=dospi/w;
  double cycle_number=6.;
  double pol_angle=6.*pi/6.;

  int nmaxt=10;//int(period*cycle_number/abs(dt));
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
  
  for(int ktime=0;ktime<=nmaxt;ktime++)
    {
      
      
      complex t = ktime*1.*dt;
      complex tt=(ktime+0.5)*1.*dt;
      
      prop_kinetic(w1, h, dt/2.);
      prop_potential(w1,h,dt);
      prop_kinetic(w1, h, dt/2.);
      w1.normalize();
	      
       /*******************************************************/
       
       fprintf(out6,"%e %e %e\n",abs(t),w1.expected_x1(),w1.expected_x2());
       cout << "\n Norm after back-forth N "<<w1.norm() <<" qN "<<w1.qnorm;
       cout << "  Fin diff <T> "<< kinetic_finite_diff(w1,h);
       cout << "  Mome space <T> "<< w1.expected_kin;
       cout << "  Potential energy "<< potential_energy(w1,h);
       cout << "  Energy "<< kinetic_finite_diff(w1,h)+potential_energy(w1,h);
      
      

	}//end time loop
    

  cout << "\n";
    cout << "  obs_kin "<< w1.expected_kin;
    //prop(w1);
cout << "  obs_kin "<< w1.expected_kin;


}//End
