#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>


class hamiltonian: public grid
{
	
public:
	
	double m1;
	double q1;		
	
	double m2;
	double q2;		
	
	double m3;
	double q3;		
	
	double a1;
	double b1;
	
	double a2;
	double b2;
	
	double a3;
	double b3;
	
	vector<double> v;  
	int something;
	
	
	/***************************************/
	// Methods.
	
	void initialize()
	{
		v.resize(n1*n2*n3, 0.);
		
		m1=1.;  //electron mass
		q1=-1;  //electron charge
		
		m2=1.;  //electron mass
		q2=-1;  //electron charge
		
		m3=1.;  //electron mass
		q3=-1;  //electron charge
		
		if(n1==1)
		{
			a1=0.;     //If the dimension is 1
			b1=0.;     //just a round zero for the operator
			
		}
		else
		{
			a1=1./2./m1;     //factor for the momentum
			b1=q1/m1;     //factor for the momentum
		}
		
		if(n2==1)
		{
			a2=0.;     //If the dimension is 1
			b2=0.;     //just a round zero for the operator
		}
		else
		{
			a2=1./2./m2;     //factor for the momentum
			b2=q2/m2;     //factor for the momentum
		}
		
		if(n3==1)
		{
			a3=0.;     //If the dimension is 1
			b3=0.;     //just a round zero for the operator
		}
		else
		{
			a3=1./2./m3;     //factor for the momentum
			b3=q3/m3;     //factor for the momentum
		}
		
	}
	
	void put_on_grid(grid g)
	{
		
		set_grid(g.n1,g.n2,g.n3,g.dx1,g.dx2,g.dx3 );
		
		symmetry_n1=g.symmetry_n1;
		symmetry_n2=g.symmetry_n2;
		symmetry_n3=g.symmetry_n3;
		
		initialize();
	}
	
	/******************************************************/
	// Fro the potential
	
	
	/******************************************************/
	// H like potential
	/******************************************************/
	void set_v_hlike3D(double charge, double soft_core)
	{
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= charge*q1/sqrt(soft_core+x1[i]*x1[i]+x2[j]*x2[j]+x3[k]*x3[k] );
				}
		
	}//End 3D pot
	
	void set_v_hlike2D(double charge, double soft_core)
	{
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					
					v[index(k,j,i)]= charge*q1/sqrt(soft_core+x1[i]*x1[i]+x2[j]*x2[j]);//+x3[k]*x3[k] );
					
				}
	}//End 2D pot
	
	void set_v_hlike1D(double charge, double soft_core)
	{
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= charge*q1/sqrt(soft_core+x1[i]*x1[i] );//+x3[k]*x3[k] );
				}  
	}//End 1D pot
	
	
	//Yukawa potential
	void set_v_hlikeY1D(double charge, double soft_core, double mp)
	{
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= charge*q1*exp(-mp*sqrt(soft_core+x1[i]*x1[i] ))/sqrt(soft_core+x1[i]*x1[i] );
				}  
	}//End 1D pot Yukawa

	
	//Yukawa potential
	void set_v_hlikeY2D(double charge, double soft_core, double mp)
	{
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= charge*q1*exp(-mp*sqrt(x2[j]*x2[j]+x1[i]*x1[i] ))/sqrt(soft_core+ x2[j]*x2[j]+ x1[i]*x1[i] );
				}  
	}//End 1D pot Yukawa	
	
	
	/******************************/
	// H like shift potential
	/******************************/
	
	void set_v_hlike_shift3D(double charge, double soft_core, double x0, double y0, double z0)
	{
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= charge*q1/sqrt(soft_core+(x1[i]-x0)*(x1[i]-x0)+(x2[j]-y0)*(x2[j]-y0)+(x3[k]-z0)*(x3[k]-z0) );
				}
		
	}//End 3D pot
	
	void set_v_hlike_shift2D(double charge, double soft_core, double x0, double y0, double z0)
	{
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= charge*q1/sqrt(soft_core+(x1[i]-x0)*(x1[i]-x0)+(x2[j]-y0)*(x2[j]-y0));//+x3[k]*x3[k] );
				}
		
	}//End 2D pot
	
	void set_v_hlike_shift1D(double charge, double soft_core, double x0, double y0, double z0)
	{
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= charge*q1/sqrt(soft_core+(x1[i]-x0)*(x1[i]-x0) );//+x3[k]*x3[k] );
				}  
	}//End 1D pot
	
	
	
	/******************************************************/
	// H2+ the potential
	/******************************************************/
	void set_hplus3D(double charge, double soft_core, double R0)
	{			
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= 
					charge*q1/sqrt(soft_core+(x1[i]-R0/2.)*(x1[i]-R0/2.)+x2[j]*x2[j]+x3[k]*x3[k] )+
					charge*q1/sqrt(soft_core+(x1[i]+R0/2.)*(x1[i]+R0/2.)+x2[j]*x2[j]+x3[k]*x3[k] );
				}
		
	}//End 3D pot		
	
	void set_hplus2D(double charge, double soft_core, double R0)
	{
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= 
					charge*q1/sqrt(soft_core+(x1[i]-R0/2.)*(x1[i]-R0/2.)+x2[j]*x2[j])+
					charge*q1/sqrt(soft_core+(x1[i]+R0/2.)*(x1[i]+R0/2.)+x2[j]*x2[j]);
				}
		
	}//End 2D pot
	
	void set_hplus1D(double charge, double soft_core, double R0)
	{					
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= 
					charge*q1/sqrt(soft_core+(x1[i]-R0/2.)*(x1[i]-R0/2.) )+
					charge*q1/sqrt(soft_core+(x1[i]+R0/2.)*(x1[i]+R0/2.) );
				}  			
	}//End 1D pot

	
	void set_hplus_shift1D(double charge, double soft_core, double R0)
	{					
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					//xshift = .8;
					v[index(k,j,i)]= 
					charge*q1/sqrt(soft_core + (x1[i]+0.)*(x1[i]+0.) )+
					charge*q1/sqrt(soft_core + ((x1[i]+0.) - R0)*((x1[i]+0.) - R0) );
				}  			
	}//End 1D pot	


	void set_co_shift_g1D(double zmol1, double zmol2, double charge, double soft_core1, double soft_core2, double R0)
	{					
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= 
					charge*zmol1*q1/sqrt(soft_core1 + x1[i]*x1[i] )+
					charge*zmol2*q1/sqrt(soft_core2 + (x1[i] - R0)*(x1[i] - R0) );
				}  			
	}//End 1D pot			

	
	
	void set_co_shift_f1D(double zmol1, double zmol2, double charge, double soft_core1, double soft_core2, double R0, double shift)
	{					
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= 
					charge*zmol1*q1/sqrt(soft_core1 + (x1[i] + shift)*(x1[i] + shift) )+
					charge*zmol2*q1/sqrt(soft_core2 + (x1[i] + shift - R0)*(x1[i] + shift - R0) );
				}  			
	}//End 1D pot				
	
	void set_co_shift1D(double zmol1, double zmol2, double charge, double soft_core, double R0)
	{					
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= 
					charge*zmol1*q1/sqrt(soft_core + x1[i]*x1[i] )+
					charge*zmol2*q1/sqrt(soft_core + (x1[i] - R0)*(x1[i] - R0) );
				}  			
	}//End 1D pot		
	

	
	//POTENTIAL WELL CO 2D
	void set_potential_CO2D(double *Zmol_i, double *Zmol_o, double *R, double *soft_core, double *sigma)
	{				
		double sqr0 = 0;
		double sqr1	 =0;
		double temp0 =0;
		double temp1 =0;		
		//cout << "\n\nZ_i1: "<< Zmol_i[0]<<"  Value From Hamiltonian Function\n\n ";
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					sqr0  = (x1[i]-R[0])*(x1[i]-R[0]) + x2[j]*x2[j];
					sqr1  = (x1[i]-R[1])*(x1[i]-R[1]) + x2[j]*x2[j];					
					
					temp0   = (Zmol_i[0] - Zmol_o[0])*exp(-sqr0/sigma[0]) + Zmol_o[0];
					temp1	= (Zmol_i[1] - Zmol_o[1])*exp(-sqr1/sigma[1]) + Zmol_o[1];					
					
					v[index(k,j,i)]= q1*( temp0/sqrt(soft_core[0]+sqr0)+
										  temp1/sqrt(soft_core[1]+sqr1));
					
				}  			
	}//End 2D pot	
	
	
	//POTENTIAL WELL CO 3D
	void set_potential_CO3D(double *Zmol_i, double *Zmol_o, double *R, double *soft_core, double *sigma)
	{				
		double sqr0 = 0;
		double sqr1	 =0;
		double temp0 =0;
		double temp1 =0;		
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					sqr0  = (x1[i]-R[0])*(x1[i]-R[0]) + x2[j]*x2[j] + x3[k]*x3[k];
					sqr1  = (x1[i]-R[1])*(x1[i]-R[1]) + x2[j]*x2[j] + x3[k]*x3[k];					
					
					temp0   = (Zmol_i[0] - Zmol_o[0])*exp(-sqr0/sigma[0]) + Zmol_o[0];
					temp1	= (Zmol_i[1] - Zmol_o[1])*exp(-sqr1/sigma[1]) + Zmol_o[1];					
					
					v[index(k,j,i)]= q1*( temp0/sqrt(soft_core[0]+sqr0)+
										 temp1/sqrt(soft_core[1]+sqr1));

				}  			
	}//End 3D pot	
	
	
	/***********************************           
	 * 	       Pot  molecular 	      *
	 ************************************/
	void moleculerN2D(double charge, double soft_core, double R0,int num_ions)
	{
		
		for(int h=0;h<num_ions;h++)
		{
			double theta=dospi*h/num_ions;
			for(int k=0;k<n3;k++)
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{ 
						v[index(k,j,i)]+= 
						charge*q1/sqrt(soft_core+
									   (x1[i]-R0*cos(theta))*(x1[i]-R0*cos(theta)) +
									   (x2[j]-R0*sin(theta))*(x2[j]-R0*sin(theta)) );
					}  
		}
		
	}//End 2D pot
	
	
	
	//POTENTIAL WELL Molecular Model for CO2 in 2D
	void set_potential_CO2_2D(double *Zmol_i, double *Zmol_o, double *R, double *soft_core, double *sigma)
	{				
		double sqr0Car		=0.;
		double sqr1Ox		=0.;
		double sqr2Ox		=0.;		
		
		double temp0Car		=0.;
		double temp1Ox		=0.;		
		double temp2Ox		=0.;				
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					sqr0Car  = (x1[i]-R[0])*(x1[i]-R[0]) + x2[j]*x2[j];
					sqr1Ox   = (x1[i]-R[1])*(x1[i]-R[1]) + x2[j]*x2[j];	
					sqr2Ox   = (x1[i]-R[2])*(x1[i]-R[2]) + x2[j]*x2[j];					
					
					
					temp0Car    = (Zmol_i[0] - Zmol_o[0])*exp( -sqr0Car/sigma[0] )  + Zmol_o[0];
					temp1Ox		= (Zmol_i[1] - Zmol_o[1])*exp( -sqr1Ox/sigma[1]  )  + Zmol_o[1];	
					temp2Ox		= (Zmol_i[2] - Zmol_o[2])*exp( -sqr2Ox/sigma[2]  )  + Zmol_o[2];		
					
					v[index(k,j,i)]= q1*( temp0Car/sqrt(soft_core[0]+sqr0Car)+
										  temp1Ox/sqrt( soft_core[1]+sqr1Ox) +
										  temp2Ox/sqrt( soft_core[2]+sqr2Ox));
					
				}  			
	}//End 2D pot		
	
	
	
	/**********************************
	 * Harmonic Oscillator             *
	 ***********************************/
	void harmonic_osc1D( double m, double w1)
	{
		
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= 
					m*w1*w1*x1[i]*x1[i]/2. ; 
					//m*w2*w2*x2[i]*x2[i]/2. +
					//m*w3*w3*x3[i]*x3[i]/2.;
				}
		
	}//End 1D pot
	
	void harmonic_osc3D( double m, double w1, double w2, double w3 )
	{			
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{ 
					v[index(k,j,i)]= 
					m*w1*w1*x1[i]*x1[i]/2. + 
					m*w2*w2*x2[i]*x2[i]/2. +
					m*w3*w3*x3[i]*x3[i]/2.;
				}			
	}//End 3D pot		
	
	/**********************************
	 *         Potential Hole          *
	 ***********************************/
	void pot_hole1D( double V0, double a0)
	{			
		for(int k=0;k<n3;k++)
			for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{  
					if (x1[i]<=-a0/2 || x1[i] >=a0/2){
						v[index(k,j,i)]= V0;
					} 
					else{
						v[index(k,j,i)]= 0;
					}  
				}			
	}//End 1D pot
	
	/**********************************
	 *         Potential Helium 1D          *
	 ***********************************/
	void set_He1D(double _charge_nuclei,double  _soft_core_ee,double _soft_core_eN )
	{			
		for(int k=0;k<n3;k++)
		    for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{  
					v[index(k,j,i)]= 
					_charge_nuclei/sqrt(x1[i]*x1[i]+_soft_core_eN)
					+_charge_nuclei/sqrt(x2[j]*x2[j]+_soft_core_eN)
					+1./sqrt((x2[j]-x1[i])*(x2[j]-x1[i])+_soft_core_ee);
					
				}			
	}//End 1D pot
	
	
	
	void save_potential(FILE *file){
		
		for(int k=0;k<n3;k++)
		    for(int j=0;j<n2;j++)
				for(int i=0;i<n1;i++)
				{  
					fprintf(file,"%e\n",v[index(k,j,i)]);
				}
	}
	
	
};

#endif
