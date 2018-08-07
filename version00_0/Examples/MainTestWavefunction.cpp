#include <iostream>
#include "grid.h"
#include "wavefunction.h"


int main()
{
  
  /****************************/
  //These are grid parameters

  int nx=20;
  int ny=10;
  int nz=300;
  
  double dx=0.1;
  double dy=0.3;
  double dz=0.2;
  /****************************/

  grid g1;                                    //Define a grid object
  g1.set_grid(nx,ny,nz,dx,dy,dz);             //Inittialize  
  
  //Sone functions of grid.

  cout << "Size in n1  "<< g1.n1 << "\n";     
  cout << "Vol elemnt  "<< g1.vol_elem() << "\n";
  cout << "Size Mb     "<< g1.aprox_size_Mb() << "\n";

  /****************************/
  // A wavefunction is a subclass of grid

  wavefunction w1;
  w1.set_grid(nx,ny,nz,dx,dy,dz);
  w1.initialize();                            //Initialize the array of the wavefunction
  cout << "Size Mb WF1   "<< w1.aprox_size_Mb() << "\n";
  
  //A wavefunction can be initialized using a grid as parameter.
  wavefunction w2;
  w2.put_on_grid(g1);
  cout << "Size Mb WF2   "<< w1.aprox_size_Mb() << "\n";

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
  // You can also create an array of functions.

  wavefunction w[10];
  for (int h=0;h<10;h++)
    {
      int nx=20+h*10;
      g1.set_grid(nx,ny,nz,dx,dy,dz);
      w[h].put_on_grid(g1);
      cout << "Size (Mb) WF"<<h<<" "<<w[h].aprox_size_Mb() << "\n";
    }

  
  
}
