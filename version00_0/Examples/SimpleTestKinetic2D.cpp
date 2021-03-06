#include <iostream>
#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "cartesian.h"

int main()
{
  
  /****************************/
  //These are grid parameters

  int nx=300;
  int ny=300;
  int nz=1;
  
  double dx=0.3;
  double dy=0.3;
  double dz=1.;
  /****************************/

  grid g1;                                    //Define a grid object
  g1.set_grid(nx,ny,nz,dx,dy,dz);             //Inittialize  
  
  //Sone info

  cout << "Size in n1  "<< g1.n1 << "\n";     
  cout << "Vol elemnt  "<< g1.vol_elem() << "\n";
  cout << "Size Mb     "<< g1.aprox_size_Mb() << "\n";

  /****************************/  
  //A wavefunction can be initialized using a grid as parameter.
  
  wavefunction w1; 
  w1.put_on_grid(g1);

  hamiltonian h;
  h.put_on_grid(g1);
  h.set_v_hlike(1., 0.6 );

  /*****************************************************/ 
  // Start the gaussian
  
   for(int k=0;k<w1.n3;k++)
    for(int j=0;j<w1.n2;j++)
      for(int i=0;i<w1.n1;i++)
	{ 	  
	  w1.w[w1.index(k,j,i)][0]= exp( -w1.x1[i]*w1.x1[i]-w1.x2[j]*w1.x2[j]-w1.x3[k]*w1.x3[k] )*cos(2.*w1.x1[i]);
	  w1.w[w1.index(k,j,i)][1]= exp( -w1.x1[i]*w1.x1[i]-w1.x2[j]*w1.x2[j]-w1.x3[k]*w1.x3[k] )*sin(2.*w1.x1[i]);;
	}

   w1.normalize(); //Normalize the wavefunction
   /*****************************************************/
   //Do the FFT in 3D

   FILE *out4,*out5,*out6,*out7;

   out4=fopen("out4.txt","w");
   out5=fopen("out5.txt","w");
   out6=fopen("out6.txt","w");
   out7=fopen("out7.txt","w");


   
   double dt=0.2;
   
   int nmaxt=200;
   int snaper1=50;//int(0.5*period/abs(dt)/2.);
   int snap=0;

   cout<<"\n nmaxt "<<nmaxt;
   
   double avect_x=0.;
   double avect_y=0.;


   int skiper1=5;
   int skiper2=5;
   int skiper3=1;
   
   for(int ktime=1;ktime<=nmaxt;ktime++)
     {
       double t = ktime*dt;
       
       prop_kinetic(w1, h, dt);
       
       /*******************************************************/
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
       
       fprintf(out6,"%e %e %e\n",t,w1.expected_x1(),w1.expected_x2());
       cout << "\n Norm after the way on and way back "<<w1.norm(); ;  
     }
     

}//End
