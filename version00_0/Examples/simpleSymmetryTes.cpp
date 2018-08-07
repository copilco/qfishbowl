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
  int nz=50.;
  
  double dx=0.4;
  double dy=0.4;
  double dz=0.4;
  /****************************/

  grid g1;                                    //Define a grid object
  g1.set_grid(nx,ny,nz,dx,dy,dz);             //Inittialize  
  
  //Sone info

  cout << "Size in n1  "<< g1.n1 << "\n";     
  cout << "Vol elemnt  "<< g1.vol_elem() << "\n";
  cout << "Size Mb     "<< g1.aprox_size_Mb() << "\n";

  /****************************/  
  //A wavefunction can be initialized using a grid as parameter.
  
  wavefunction w1,w2,w3,w4; 
  w1.put_on_grid(g1);
  w2.put_on_grid(g1);
  w3.put_on_grid(g1);
  w4.put_on_grid(g1);

  hamiltonian h;
  h.put_on_grid(g1);

  double charge_nuclei=1.;
  double soft_core=0.5;//2.;//0.6;

  h.set_v_hlike(charge_nuclei, soft_core );

  /*****************************************************/ 
  // Start the gaussian
  gaussian(w1);
  gaussian(w2);
  gaussian(w3);
  gaussian(w4);
  
  w1.normalize(); //Normalize the wavefunction
  w2.normalize(); //Normalize the wavefunction
  w3.normalize(); //Normalize the wavefunction
  w4.normalize(); //Normalize the wavefunction
  
  /*****************************************************/
  //Do the FFT in 3D
  
  FILE *out0,*out1,*out2,*out3,*out4,*out5,*out6,*out7;
  
  out0=fopen("3exout0.txt","w");
  out1=fopen("3exout1.txt","w");
  out2=fopen("3exout2.txt","w");
  out3=fopen("3exout3.txt","w");  
  out4=fopen("3exout4.txt","w");
  out5=fopen("3exout5.txt","w");
  out6=fopen("3exout6.txt","w");
  out7=fopen("3exout7.txt","w");

 
 /*************************************************************************************************/
  complex dt=0.05;//-I*0.05;
  
  
  double e0=sqrt(1.e14/3.5e16);
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
  int skiper1=1;
  int skiper2=1;
  int skiper3=1;
  
  dt=-I*0.05;
  nmaxt=5000;
  double e2=1.;


  /*******************************************************************************************/
  for(int ktime=0;ktime<=nmaxt;ktime++)
    {
      
      complex t = ktime*1.*dt;
      complex tt=(ktime+0.5)*1.*dt;
      
      prop_kinetic(w1, h, dt/2.);
      prop_potential(w1,h,dt);
      prop_kinetic(w1, h, dt/2.);
    

    
      
      w1.normalize();

      /*******************************************************/
      //Observables      
      if((ktime%500)==0)
	{
      	  cout << "\n N0 "<<w1.norm() <<" "<<w1.q_norm();
	  cout << "  E "<< kinetic_finite_diff(w1,h)+potential_energy(w1,h);
	  
	  double e1=kinetic_finite_diff(w1,h)+potential_energy(w1,h);
	  double err=log10(abs(e1-e2));
	  cout << "  err " << err;
	  e2=e1;
	  
	}
      /*******************************************************/

    }
  
  for(int k=0;k<w1.n3/skiper3;k++)
    for(int j=0;j<w1.n2/skiper2;j++)
      for(int i=0;i<w1.n1/skiper1;i++)
	{
	  double norm=
	    w1.dx1*w1.dx2*w1.dx3*(w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				  w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]+
				  w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				  w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]);
	  fprintf(out5,"%e\n",norm);
	}


  /*******************************************************************************************/
  



  /*******************************************************************************************/
  int nmaxt2=5*nmaxt;
  for(int ktime=0;ktime<=nmaxt2;ktime++)
    {
      
      complex t = ktime*1.*dt;
      complex tt=(ktime+0.5)*1.*dt;
      
      prop_kinetic(w2, h, dt/2.);
      prop_potential(w2,h,dt);
      prop_kinetic(w2, h, dt/2.);
      
      w2.normalize();
      project_out(w2,w1);
      w2.normalize();

      /*
      ///This symmetry constrain makes the wavefunction 
      // to be antisymmetric with respect to the diagonal
      // Creates the 2P state without haveing the @s built.
      for(int k=0;k<w2.n3;k++)
	for(int j=0;j<w2.n2;j++)
	    {
	      w2.w[w2.index(k,j,j)][0]=0;
	      w2.w[w2.index(k,j,j)][1]=0;
	    }
      
      /*******************************************************/
      //Observables
      
      if((ktime%50)==0)
	{
	  cout << "\n N1 "<<w2.norm() <<" "<<w2.q_norm();
	  cout << "  E "<< kinetic_finite_diff(w2,h)+potential_energy(w2,h);
	  cout << "  Projection " << norm(projection(w2,w1));
	  
	  double e1=kinetic_finite_diff(w2,h)+potential_energy(w2,h);
	  double err=log10(abs(e1-e2));
	  
	  if(err<-10.)
	    ktime=nmaxt2;
	  
	  cout << "  err " << err;
	  e2=e1;
	}
      /*******************************************************/
      
    }//End of time loop
  
  
  for(int k=0;k<w2.n3/skiper3;k++)
    for(int j=0;j<w2.n2/skiper2;j++)
      for(int i=0;i<w2.n1/skiper1;i++)
	{
	  
	  double norm=
	    w2.dx1*w2.dx2*w2.dx3*(w2.w[w2.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				  w2.w[w2.index(k*skiper3,j*skiper2,i*skiper1)][0]+
				  w2.w[w2.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				  w2.w[w2.index(k*skiper3,j*skiper2,i*skiper1)][1]);
	  fprintf(out5,"%e\n",norm);
	}


  /*******************************************************************************************/

}//End
