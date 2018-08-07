#include <iostream>
#include "grid.h"

int main()
{
  
  int nx=100;
  int ny=1000;
  int nz=300;
  
  double dx=0.1;
  double dy=0.3;
  double dz=0.2;

  grid g1;
  
  g1.set_grid(nx,ny,nz,dx,dy,dz);
  
  cout << "Size in n1  "<< g1.n1 << "\n";
  cout << "Vol elemnt  "<< g1.vol_elem() << "\n";
  cout << "Size Mb     "<< g1.aprox_size_Mb() << "\n";
  
}
