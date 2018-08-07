#include <iostream>
#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "cartesian.h"
using namespace std;

int main()
{
	
	/****************************/
	//These are grid parameters
	
	int nx=800;
	int ny=1;
	int nz=1.;
	
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
	
	double charge_nuclei=1.;
	double soft_core=2.;
	//h.set_v_hlike1D(charge_nuclei, soft_core );
	//h.set_v_hlike_shift1D(charge_nuclei, soft_core, 10.,0,0);
	h.set_hplus1D(charge_nuclei, soft_core, 12.);
	
	/*****************************************************/ 
	// Start the gaussian
	gaussian(w1);  
	
	w1.normalize(); //Normalize the wavefunction
	
	/*****************************************************/
	//Files
	
	FILE *out0,*out1,*out2,*out3,*out4,*out5,*out6,*out7;
	
	out0=fopen("out0.txt","w");
	out1=fopen("out1.txt","w");
	out2=fopen("out2.txt","w");
	out3=fopen("out3.txt","w");  
	out4=fopen("out4.txt","w");
	out5=fopen("out5.txt","w");
	out6=fopen("out6.txt","w");
	out7=fopen("out7.txt","w");
	
	for (int i=0;i<g1.n1;i++)
		fprintf(out0,"%e\n",g1.x1[i]);

	
	/*************************************************************************************************/
	complex dt=-I*0.05;

	
	/*************************************************************************************************/
	
	for(int ktime=0;ktime<=10000;ktime++)
    {

		prop_kinetic(w1, h, dt/2.);
		prop_potential(w1,h,dt);
		prop_kinetic(w1, h, dt/2.);
		w1.normalize();
		
		if((ktime%200)==0)
			cout << "\nImaginary loop converging Energy .."<< kinetic_finite_diff(w1,h)+potential_energy(w1,h);
		
		if((ktime%4)==0)
		{
			fprintf(out4,"%e\n",ktime*abs(dt));
			snapshot(out5,w1,1,1,1);
		}
		
    }//End imaginary loop
	
	
		
}//End
