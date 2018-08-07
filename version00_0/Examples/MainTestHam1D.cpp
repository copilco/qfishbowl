#include <iostream>
#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"


int main()
{
  
  /****************************/
  //These are grid parameters

  int nx=400;
  int ny=1;
  int nz=1;
  
  double dx=0.3;
  double dy=1.;
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
	  w1.w[w1.index(k,j,i)][0]= exp( -w1.x1[i]*w1.x1[i]-w1.x2[j]*w1.x2[j]-w1.x3[k]*w1.x3[k] )*cos(w1.x1[i]);
	  w1.w[w1.index(k,j,i)][1]= exp( -w1.x1[i]*w1.x1[i]-w1.x2[j]*w1.x2[j]-w1.x3[k]*w1.x3[k] )*sin(w1.x1[i]);;
	}

   w1.normalize(); //Normalize the wavefunction
   /*****************************************************/
   //Do the FFT in 3D

   FILE *out4,*out5,*out6;

   out4=fopen("out4.txt","w");
   out5=fopen("out5.txt","w");
   out6=fopen("out6.txt","w");


   double factor=1./w1.n1/w1.n2/w1.n3;     
   double kinetic;
   double dt=0.2;

   int skiper1=1;
   int skiper2=1;
   int skiper3=1;
   
   for(int ktime=1;ktime<=200;ktime++)
     {

       
       /**********************/
       fftw_execute(w1.p3DF);
       /**********************/
              
       for(int k=0;k<w1.n3;k++)
	 for(int j=0;j<w1.n2;j++)
	   for(int i=0;i<w1.n1;i++)
	     {
	       
	       kinetic=
		 w1.q1[i]*w1.q1[i]*h.a1+
		 w1.q2[j]*w1.q2[j]*h.a2+
		 w1.q3[k]*w1.q3[k]*h.a3;
	       
	       //expected_kin+=dqz*dqx*dqy*kinetic*
	       //  (phic[index(k,j,i)][0]*phic[index(k,j,i)][0]+phic[index(k,j,i)][1]*phic[index(k,j,i)][1])
	       //  *(dx*dy*dz*dx*dy*dz/dospi/dospi/dospi);
	       
	       complex kinoperator=exp(-I*dt*kinetic );
	       
	       double aux0r=w1.w[w1.index(k,j,i)][0];
	       double aux0i=w1.w[w1.index(k,j,i)][1];
	       
	       w1.w[w1.index(k,j,i)][0]=( aux0r*real(kinoperator)+aux0i*imag(kinoperator) )*factor;
	       w1.w[w1.index(k,j,i)][1]=( aux0i*real(kinoperator)-aux0r*imag(kinoperator) )*factor;
	       	       
	     }
       

       /**********************/
       fftw_execute(w1.p3DB);
       /**********************/

       for(int k=0;k<w1.n3/skiper3;k++)
	 for(int j=0;j<w1.n2/skiper2;j++)
	   for(int i=0;i<w1.n1/skiper1;i++)
	     {
	       
	       double norm=
		 w1.dx1*w1.dx2*w1.dx3*(w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]+
				       w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]);
	       fprintf(out5,"%e\n",norm);
	     }


       
       fprintf(out6,"%d %e\n",ktime,w1.expected_x1());
       cout << "\n Norm after the way on and way back "<<w1.norm() << " norm with factor " << w1.norm()*factor*factor;  
       //w1.normalize();
     }
   

}
