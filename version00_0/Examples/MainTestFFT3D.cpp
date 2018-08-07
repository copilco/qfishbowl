#include <iostream>
#include "grid.h"
#include "wavefunction.h"


int main()
{
  
  /****************************/
  //These are grid parameters

  int nx=40;
  int ny=40;
  int nz=40;
  
  double dx=0.3;
  double dy=0.3;
  double dz=0.3;
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

  cout << "Size in n1  "<< w1.n1 << "\n";     
  cout << "Vol elemnt  "<< w1.vol_elem() << "\n";
  cout << "Size Mb     "<< w1.aprox_size_Mb() << "\n";

  //Write to file the axes in coordinate and momentum space
  FILE *out1,*out2,*out3;
  out1 = fopen("axis1.txt","w");
  out2 = fopen("axis2.txt","w");
  out3 = fopen("axis3.txt","w");
  
  cout << "Space spanned in q1 (a.u.) WF1   "<< w1.sizeq1_au() << "\n";
  cout << "Space spanned in x1 (a.u.) WF1   "<< w1.size1_au() << "\n";

  cout << "Space spanned in q2 (a.u.) WF1   "<< w1.sizeq2_au() << "\n";
  cout << "Space spanned in x2 (a.u.) WF1   "<< w1.size2_au() << "\n";

  cout << "Space spanned in q3 (a.u.) WF1   "<< w1.sizeq3_au() << "\n";
  cout << "Space spanned in x3 (a.u.) WF1   "<< w1.size3_au() << "\n";

  for(int h=0;h<w1.n1;h++)
    fprintf(out1,"%e %e\n",w1.x1[h],w1.q1[h]);

  for(int h=0;h<w1.n2;h++)
    fprintf(out2,"%e %e\n",w1.x2[h],w1.q2[h]);

  for(int h=0;h<w1.n3;h++)
    fprintf(out3,"%e %e\n",w1.x3[h],w1.q3[h]);

  /*****************************************************/ 
  // Start the gaussian

   for(int k=0;k<w1.n3;k++)
    for(int j=0;j<w1.n2;j++)
      for(int i=0;i<w1.n1;i++)
	{ 	  
	  w1.w[w1.index(k,j,i)][0]= exp( -w1.x1[i]*w1.x1[i]-w1.x2[j]*w1.x2[j]-w1.x3[k]*w1.x3[k] );
	  w1.w[w1.index(k,j,i)][1]=0.;
	}

   w1.normalize(); //Normalize the wavefunction
   /*****************************************************/ 
   //Do the FFT in 3D

   FILE *out4,*out5,*out6;

   out4=fopen("fileinX.txt","w");
   out5=fopen("fileinQ.txt","w");
   out6=fopen("fileoutX.txt","w");

   int skiper=1;
   for(int k=0;k<w1.n3;k++)
     for(int j=0;j<w1.n2;j++)
       for(int i=0;i<w1.n1;i++)
	 {
	   double norm=
	     w1.dx1*w1.dx2*w1.dx3*(w1.w[w1.index(k*skiper,j*skiper,i*skiper)][0]*w1.w[w1.index(k*skiper,j*skiper,i*skiper)][0]
			  +w1.w[w1.index(k*skiper,j*skiper,i*skiper)][1]*w1.w[w1.index(k*skiper,j*skiper,i*skiper)][1]);
	   fprintf(out4,"%e\n",norm);
	 }


   double factor=1./w1.n1/w1.n2/w1.n3;  
   /**********************/
   fftw_execute(w1.p3DF);
   /**********************/
   
   for(int k=0;k<w1.n3;k++)
     for(int j=0;j<w1.n2;j++)
       for(int i=0;i<w1.n1;i++)
	 {
	   double norm=
	     w1.dq1*w1.dq2*w1.dq3*(w1.w[w1.index(k*skiper,j*skiper,i*skiper)][0]*w1.w[w1.index(k*skiper,j*skiper,i*skiper)][0]
			     +w1.w[w1.index(k*skiper,j*skiper,i*skiper)][1]*w1.w[w1.index(k*skiper,j*skiper,i*skiper)][1]);
	   fprintf(out5,"%e\n",norm);
	 }
   
   /**********************/
   fftw_execute(w1.p3DB);
   /**********************/
   
   
   for(int k=0;k<w1.n3;k++)
     for(int j=0;j<w1.n2;j++)
       for(int i=0;i<w1.n1;i++)
	 {
	   double norm=
	   factor*factor* w1.dx1*w1.dx2*w1.dx3*(w1.w[w1.index(k*skiper,j*skiper,i*skiper)][0]*w1.w[w1.index(k*skiper,j*skiper,i*skiper)][0]
			  +w1.w[w1.index(k*skiper,j*skiper,i*skiper)][1]*w1.w[w1.index(k*skiper,j*skiper,i*skiper)][1]);
	   fprintf(out6,"%e\n",norm);
	 }

   cout << "\n Norm after the way on and way back "<<w1.norm() << " norm with factor " << w1.norm()*factor*factor;  

}
