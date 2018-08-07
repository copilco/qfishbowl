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
							
		/******************************************************/
		// H like shift potential
		/******************************************************/
		
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
		
		/******************************
		* Harmonic Oscillator
		******************************/
		void harmonic_osc1D(double charge, double m, double w0)
		{
			
			for(int k=0;k<n3;k++)
				for(int j=0;j<n2;j++)
					for(int i=0;i<n1;i++)
					{ 
						v[index(k,j,i)]= 
						m*w0*w0*x1[i]*x1[i]/2.;
					}
			
		}//End 1D pot
		
		
	};
