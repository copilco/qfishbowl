/* 
 * Set of Bound States with the Split Operator Tecnique on imaginary propagation time
 */

#include <iostream>
#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "cartesian.h"
#include "constant.h"
#include<sstream>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

int main()
{
	/****************************
	 *   PARAMETERS ATOM 2D      *
	 ****************************/
	//      mkdir("BoundStateSOHe",0777);
	int Nw   = 4;               //Number of Bound States
	int Nout = 2*(Nw+1);
	FILE *bin[Nw];
	FILE *s[Nw];	
	FILE *out[Nout];
	
	string abin="bstate";
	string aste="state";
	string aout="out";
	for(int i=0;i<Nout;i++)
	{
		if(i<Nw){
			stringstream bbin;
			stringstream bste;
			bbin << abin << i;
			bste << aste << i;
			string afile = bbin.str()+".bin";
			string bfile = bste.str()+".txt";  
			bin[i] = fopen(afile.c_str(),"wb");
			s[i] = fopen(bfile.c_str(),"w");
		}
		stringstream bout;
		bout << aout << i;
		string cfile = bout.str()+".txt";
		out[i] = fopen(cfile.c_str(),"w");
	}
	
	//Parameter of grid
	int nx    = 4000;//6000;
	int ny    = 1;
	int nz    = 1;
    
	double dx = 0.05;
	double dy = 1.;
	double dz = 1.; 	
	
	
	
	//Potential parameters and imaginary loop
	int Ntime   = 15000;     		     //Iteration loop imaginary propagation Optimum 1620
	int snap2   = 4;
	int snap3   = 4;
	
	
	
	//Complex time for imaginary time evolution 
	complex dt  = 0.008;			     //complex time	
	
	
	
	//Parameter of atom Helio+
	double charge_nuclei = 1.;		       //Charge nuclei
	double soft_core     = 2.;	           //Soft_core parameter soft_core=0.487688
	
	
	
	
	//snapers
	int snaper1  = int(0.125*110);
	int snaper2  = 100;
	
	
	
	//Objet grid
	grid g1;
	g1.set_grid(nx,ny,nz,dx,dy,dz);  //Initialize grid
	XQsnapshot(out[0], g1, 1, 1, 1);  
	
	
	
	
	//A wavefunction can be initialized using a grid as parameter.    
	wavefunction w[Nw];	
	for(int j=0;j<Nw;j++)     
		w[j].put_on_grid(g1);
	
	
	
	
	//Start the gaussian
	for (int j=0;j<Nw;j++)
	{
		gaussian(w[j]);
		w[j].normalize();
	}
	
	
	
	
	
	//Start hamiltonian's potential grid "minor"
	hamiltonian h1;
	h1.put_on_grid(g1);
	h1.set_v_hlike1D(charge_nuclei, soft_core); //atom potential    
	
	
	
	
	
	
	/****************************
	 *     STATES CALCULATED	  *
	 ****************************/
	dt=-I*dt;
	double TcalValue =-0.5; 	       //energy of state 1s
	double Energ0    = TcalValue;	   //energy of state 1s
	double Energ1    = 0; 
	double Energ2    = 0;		
	double Energ3    = 0;
	double Energ4    = 0;        
	int kont         = 0;
	
	
	
	
	for (int jw=0;jw<Nw;jw++)
	{		
		
		
		
		//Start loop for calculation the jw-th state 
		for(int ktime=0;ktime<=Ntime;ktime++)
		{   
			
			
			prop_kinetic(w[jw], h1, dt/2.);
			prop_potential(w[jw], h1, dt);
			prop_kinetic(w[jw], h1, dt/2.);
			w[jw].normalize();
			
			
			
			for(int jp=0;jp<jw;jp++)
			{                     
				project_out(w[jw],w[jp]);                         
				w[jw].normalize();
			}		
			
			
			
			
			if((ktime%(Ntime/15))==0)
			{
				double valfinite=kinetic_finite_diff(w[jw],h1)+potential_energy(w[jw],h1);
				Energ4=q_expected_kinetic(w[jw],h1)+potential_energy(w[jw],h1);
				double ThError=abs( Energ4-Energ0 );
				double ConvError=abs( Energ4-TcalValue );
				
				cout << "\n N"<< jw <<":  "<<w[jw].norm();
				cout << "   E"<< jw <<":  "<<valfinite;
				cout << "   E"<< jw <<":  "<<Energ4;
				cout << "   ThError"<< jw <<":  "<< ThError;
				cout << "   ConvError"<< jw <<":  "<< ConvError;
				fprintf(out[2*jw+1],"%e %10.20e\n", ConvError, w[jw].expected_absx1());
				snapshot(out[2*jw+2],w[jw],1,1,1);
				if (abs(ConvError)<=1.0e-15){
					kont=kont+1;
					if (kont==7) break;
				}
			}
			TcalValue=q_expected_kinetic(w[jw],h1)+potential_energy(w[jw],h1);
			for(int j=0;j<w[jw].n1*w[jw].n2; j++)        
				w[jw].w[j][1]=0.0;
		}//end loop state jw
		
		
		
		fprintf(s[jw],"%e %e \n", q_expected_kinetic(w[jw],h1)+potential_energy(w[jw],h1), 0.000);
		
		for (int j=0;j<w[jw].n1*w[jw].n2*w[jw].n3;j++)
			fprintf(s[jw],"%e %e \n", w[jw].w[j][0], w[jw].w[j][1]);
		
		
		biwrite(bin[jw], w[jw]);
		Ntime+=floor(Ntime/2)+1;
		kont=0;
	}//Set State
	
	
	
	
}//End
