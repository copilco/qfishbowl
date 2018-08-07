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
	h.set_v_hlike1D(charge_nuclei, soft_core );
	
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
	
	double e0=sqrt(1.e14/3.5e16);
	double w=0.057;
	double period=dospi/w;
	double cycle_number=6.;
	complex dt=-I*period/5000.;
	
	int nmaxt=int(period*cycle_number/abs(dt));
	int snaper1=int(0.125*period/abs(dt));
	int snaper2=1;//int(0.125*period/abs(dt));
	
	/*************************************************************************************************/
	
	for(int ktime=0;ktime<=10000;ktime++)
    {

		prop_kinetic(w1, h, dt/2.);
		prop_potential(w1,h,dt);
		prop_kinetic(w1, h, dt/2.);
		w1.normalize();
		
		if((ktime%100)==0)
			cout << "\nImaginary loop converging Energy .."<< kinetic_finite_diff(w1,h)+potential_energy(w1,h);
		
    }//End imaginary loop
	
	
	//Real time
	int skiper1=5;
	int skiper2=1;
	int skiper3=1;
	
	dt=abs(dt);
	double avect_x=0.;
	double avect_y=0.;
	
	nmaxt=int(period*cycle_number/abs(dt));
	for(int ktime=0;ktime<=nmaxt;ktime++)
    {
		
		complex t = ktime*1.*dt;
		complex tt=(ktime+0.5)*1.*dt;
		
		
		double efield_x=e0*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t)+pi/2.);
		
		avect_x+=-efield_x*abs(dt)*lightC_au;
		
		fprintf(out3,"%e %e %e\n",abs(t),efield_x,avect_x/lightC_au);
		
		prop_kinetic_laser_AP(w1, h, dt/2., avect_x,0.,0.);
		prop_potential(w1,h,dt);
		prop_kinetic_laser_AP(w1, h, dt/2., avect_x,0.,0.);
		
		
		/*******************************************************/
		if((ktime%4)==0)
		{
			fprintf(out4,"%e\n",abs(t));
			snapshot(out5,w1,1,1,1);
		}
		

		/*******************************************************/
		//Observables

		fprintf(out6,"%e %e\n",abs(t),w1.norm());
		
		if((ktime%100)==0)
			cout << "\n Norm "<<w1.norm() <<" "<<w1.q_norm() <<"  Energy "<< kinetic_finite_diff(w1,h)+potential_energy(w1,h);
		
		double dipole_para=0.;
		double dipole_perp=0.;
		
		double eldipole_para=0.;
		double eldipole_perp=0.;
		
		double dipole_x_axis=0.;
		double dipole_y_axis=0.;
		
		skiper1=skiper2=skiper3=1;
		
		for(int k=0;k<w1.n3/skiper3;k++)
			for(int j=0;j<w1.n2/skiper2;j++)
				for(int i=0;i<w1.n1/skiper1;i++)
				{		     
					double dipole_x=charge_nuclei*pow(soft_core+(h.x1[i])*(h.x1[i]),-3./2.)*h.x1[i];
					double norm_temp=(w1.w[w1.index(k,j,i)][0]*
									  w1.w[w1.index(k,j,i)][0]+
									  w1.w[w1.index(k,j,i)][1]*
									  w1.w[w1.index(k,j,i)][1]);
					dipole_x_axis+=w1.dx1*w1.dx2*w1.dx3*(dipole_x)*norm_temp;
					
				}
		
		fprintf(out7,"%e %e \n",abs(t),dipole_x_axis);
		
		
		absorber(w1,0.1,0.,0.,0.1,0.,0.,1./6.);
		
    }//end time loop
	
}//End
