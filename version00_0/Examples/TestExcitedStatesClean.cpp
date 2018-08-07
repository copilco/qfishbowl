//#include <iostream>
#include <iostream>
#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "cartesian.h"
//#include "constant.h"
using namespace std;
//#include <complex.h>
//#define complex complex<double>

int main()
{
  
  /****************************/
  //These are grid parameters

  int nx=100;
  int ny=100;
  int nz=1.;
  
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
  
  wavefunction w1,w2,w3,w4; 
  w1.put_on_grid(g1);
  w2.put_on_grid(g1);
  w3.put_on_grid(g1);
  w4.put_on_grid(g1);

  hamiltonian h;
  h.put_on_grid(g1);

  double charge_nuclei=1.;
  double soft_core=0.5;//2.;//0.6;

  h.set_v_hlike(charge_nuclei, soft_core );

  /*****************************************************/ 
  // Start the gaussian
  gaussian(w1);
  gaussian(w2);
  gaussian(w3);
  gaussian(w4);
  
  w1.normalize(); //Normalize the wavefunction
  w2.normalize(); //Normalize the wavefunction
  w3.normalize(); //Normalize the wavefunction
  w4.normalize(); //Normalize the wavefunction
  
  /*****************************************************/
  //Do the FFT in 3D
  
  FILE *out0,*out1,*out2,*out3,*out4,*out5,*out6,*out7;
  
  out0=fopen("exout0.txt","w");
  out1=fopen("exout1.txt","w");
  out2=fopen("exout2.txt","w");
  out3=fopen("exout3.txt","w");  
  out4=fopen("exout4.txt","w");
  out5=fopen("exout5.txt","w");
  out6=fopen("exout6.txt","w");
  out7=fopen("exout7.txt","w");

 
 /*************************************************************************************************/
  complex dt=0.05;//-I*0.05;
  
  
  double e0=sqrt(1.e14/3.5e16);
  double w=0.057;
  double period=dospi/w;
  double cycle_number=6.;
  double pol_angle=6.*pi/6.;

  int nmaxt=int(period*cycle_number/abs(dt));
  int snaper1=int(0.125*period/abs(dt));

  cout<<"\n nmaxt "<<nmaxt;
  
  double avect_x=0.;
  double avect_y=0.;

  int snap=0;
  cout<<"\n nmaxt "<<nmaxt;
  /*************************************************************************************************/
  int skiper1=1;
  int skiper2=1;
  int skiper3=1;
  
  dt=-I*0.05;
  nmaxt=3000;
  double e2=1.;
  /*******************************************************************************************/
  for(int ktime=0;ktime<=nmaxt;ktime++)
    {
      
      complex t = ktime*1.*dt;
      complex tt=(ktime+0.5)*1.*dt;
      
      prop_kinetic(w1, h, dt/2.);
      prop_potential(w1,h,dt);
      prop_kinetic(w1, h, dt/2.);
      w1.normalize();

      /*******************************************************/
      //Observables      
      if((ktime%500)==0)
	{
      	  cout << "\n N0 "<<w1.norm() <<" "<<w1.q_norm();
	  cout << "  E "<< kinetic_finite_diff(w1,h)+potential_energy(w1,h);
	  
	  double e1=kinetic_finite_diff(w1,h)+potential_energy(w1,h);
	  double err=log10(abs(e1-e2));
	  cout << "  err " << err;
	  e2=e1;
	  
	}
      /*******************************************************/

    }
  
  for(int k=0;k<w1.n3/skiper3;k++)
    for(int j=0;j<w1.n2/skiper2;j++)
      for(int i=0;i<w1.n1/skiper1;i++)
	{
	  double norm=
	    w1.dx1*w1.dx2*w1.dx3*(w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				  w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]+
				  w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				  w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]);
	  fprintf(out5,"%e\n",norm);
	}


  /*******************************************************************************************/
  int nmaxt2=5*nmaxt;
  for(int ktime=0;ktime<=nmaxt2;ktime++)
    {
      
      complex t = ktime*1.*dt;
      complex tt=(ktime+0.5)*1.*dt;
      
      prop_kinetic(w2, h, dt/2.);
      prop_potential(w2,h,dt);
      prop_kinetic(w2, h, dt/2.);
      
      w2.normalize();
      project_out(w2,w1);
      w2.normalize();

      /*******************************************************/
      //Observables
      
      if((ktime%50)==0)
	{
	  cout << "\n N1 "<<w2.norm() <<" "<<w2.q_norm();
	  cout << "  E "<< kinetic_finite_diff(w2,h)+potential_energy(w2,h);
	  cout << "  Projection " << norm(projection(w2,w1));
	  
	  double e1=kinetic_finite_diff(w2,h)+potential_energy(w2,h);
	  double err=log10(abs(e1-e2));

	  if(err<-10.)
	    ktime=nmaxt2;

	  cout << "  err " << err;
	  e2=e1;
	}
      /*******************************************************/

    }//End of time loop


  for(int k=0;k<w2.n3/skiper3;k++)
    for(int j=0;j<w2.n2/skiper2;j++)
      for(int i=0;i<w2.n1/skiper1;i++)
	{
	  
	  double norm=
	    w2.dx1*w2.dx2*w2.dx3*(w2.w[w2.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				  w2.w[w2.index(k*skiper3,j*skiper2,i*skiper1)][0]+
				  w2.w[w2.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				  w2.w[w2.index(k*skiper3,j*skiper2,i*skiper1)][1]);
	  fprintf(out5,"%e\n",norm);
	}


  /*******************************************************************************************/
  int nmaxt3=5*nmaxt;
  for(int ktime=0;ktime<=nmaxt3;ktime++)
    {
      
      complex t = ktime*1.*dt;
      complex tt=(ktime+0.5)*1.*dt;
      
      prop_kinetic(w3, h, dt/2.);
      prop_potential(w3,h,dt);
      prop_kinetic(w3, h, dt/2.);
      
      /*
      for(int k=0;k<w2.n3;k++)
	for(int j=0;j<w2.n2/2;j++)
	  for(int i=0;i<w2.n1;i++)
	    {
	      w2.w[w2.index(k,j,i)][0]=-(w2.w[w2.index(k,(w2.n2-1)-j,i)][0]);
	      w2.w[w2.index(k,j,i)][1]=-(w2.w[w2.index(k,(w2.n2-1)-j,i)][1]);
	    } 
      */

      w3.normalize();
      project_out(w3,w1);
      w3.normalize();
      project_out(w3,w2);
      w3.normalize();

      /*******************************************************/
      //Observables
      
      if((ktime%50)==0)
	{
	  
    
	  cout << "\n N2 "<<w3.norm() <<" "<<w3.q_norm();
	  cout << "  E "<< kinetic_finite_diff(w3,h)+potential_energy(w3,h);
	  cout << "  Projection " << norm(projection(w3,w2));
	  double e1=kinetic_finite_diff(w3,h)+potential_energy(w3,h);
	  double err=log10(abs(e1-e2));
	  
	  if(err<-14.)
	    ktime=nmaxt3;
	  
	  cout << "  err " << err;
	  e2=e1;
	}
      /*******************************************************/

    }


  for(int k=0;k<w3.n3/skiper3;k++)
    for(int j=0;j<w3.n2/skiper2;j++)
      for(int i=0;i<w3.n1/skiper1;i++)
	{
	  
	  double norm=
	    w3.dx1*w3.dx2*w3.dx3*(w3.w[w3.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				  w3.w[w3.index(k*skiper3,j*skiper2,i*skiper1)][0]+
				  w3.w[w3.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				  w3.w[w3.index(k*skiper3,j*skiper2,i*skiper1)][1]);
	  fprintf(out5,"%e\n",norm);
	}


  /*******************************************************************************************/
  int nmaxt4=5*nmaxt;
  for(int ktime=0;ktime<=nmaxt4;ktime++)
    {
      
      complex t = ktime*1.*dt;
      complex tt=(ktime+0.5)*1.*dt;
      
      prop_kinetic(w4, h, dt/2.);
      prop_potential(w4,h,dt);
      prop_kinetic(w4, h, dt/2.);
      
      /*
      for(int k=0;k<w2.n3;k++)
	for(int j=0;j<w2.n2/2;j++)
	  for(int i=0;i<w2.n1;i++)
	    {
	      w2.w[w2.index(k,j,i)][0]=-(w2.w[w2.index(k,(w2.n2-1)-j,i)][0]);
	      w2.w[w2.index(k,j,i)][1]=-(w2.w[w2.index(k,(w2.n2-1)-j,i)][1]);
	    } 
      */

      w4.normalize();
      project_out(w4,w1);
      w4.normalize();
      project_out(w4,w2);
      w4.normalize();
      project_out(w4,w3);
      w4.normalize();

      /*******************************************************/
      //Observables
      
      if((ktime%50)==0)
	{
	  
	  
	  cout << "\n N3 "<<w4.norm() <<" "<<w4.q_norm();
	  cout << "  E "<< kinetic_finite_diff(w4,h)+potential_energy(w4,h);
	  cout << "  Projection " << norm(projection(w4,w3));

	  double e1=kinetic_finite_diff(w4,h)+potential_energy(w4,h);
	  double err=log10(abs(e1-e2));
	  
	  if(err<-14.)
	    ktime=nmaxt4;
	  
	  cout << "  err " << err;
	  e2=e1;

	}
      /*******************************************************/
    }//end time loop


  for(int k=0;k<w4.n3/skiper3;k++)
    for(int j=0;j<w4.n2/skiper2;j++)
      for(int i=0;i<w4.n1/skiper1;i++)
	{
	  
	  double norm=
	    w4.dx1*w4.dx2*w4.dx3*(w4.w[w4.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				  w4.w[w4.index(k*skiper3,j*skiper2,i*skiper1)][0]+
				  w4.w[w4.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				  w4.w[w4.index(k*skiper3,j*skiper2,i*skiper1)][1]);
	  fprintf(out5,"%e\n",norm);
	}




  /*******************************************************************************************/

  //Real time
  dt=0.05;
  nmaxt=int(period*cycle_number/abs(dt));
  for(int ktime=0;ktime<=-nmaxt;ktime++)
    {
      
      complex t = ktime*1.*dt;
      complex tt=(ktime+0.5)*1.*dt;

       
      double efield_x=e0*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t)+pi/2.);
      double efield_y=e0*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t)/2./cycle_number)*sin(w*abs(t)+pi/2.);
      
      avect_x+=-efield_x*abs(dt)*lightC_au;
      avect_y+=-efield_y*abs(dt)*lightC_au;
      
      fprintf(out3,"%e %e %e %e %e\n",abs(t),efield_x,avect_x/lightC_au,efield_y,avect_y/lightC_au);
	 
      prop_kinetic_laser_AP(w1, h, dt, avect_x,0.,0.);
      //prop_kinetic(w1, h, dt/2.);
      prop_potential(w1,h,dt);
      prop_kinetic_laser_AP(w1, h, dt, avect_x,0.,0.);
      //prop_kinetic(w1, h, dt/2.);

      /*
      prop_kinetic_laser_AP(w1, h, dt, avect_x,0.,0.);
      //prop_kinetic(w2, h, dt/2.);
      prop_potential(w2,h,dt);
      prop_kinetic_laser_AP(w1, h, dt, avect_x,0.,0.);
      //prop_kinetic(w2, h, dt/2.);
      */

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
      //Observables
      
      fprintf(out6,"%e %e %e\n",abs(t),w1.expected_x1(),w1.expected_x2());
      cout << "\n N "<<w1.norm() <<" "<<w1.q_norm();
      cout << "  K E "<< kinetic_finite_diff(w1,h);
      cout << "  P E "<< potential_energy(w1,h);
      cout << "  E "<< kinetic_finite_diff(w1,h)+potential_energy(w1,h);
      cout << "  kin "<< q_expected_kinetic(w1,h);
      //cout << "  dq "<< w1.dq1 << " "<< w1.dq2 << " "<< w1.dq3;
	  //cout << "  f "<< w1.kfactor1 << " "<< w1.kfactor2 << " "<< w1.kfactor3;
	 
      cout << " *N "<<w2.norm() <<" "<<w2.q_norm();
      cout << "  K E "<< kinetic_finite_diff(w2,h);
      cout << "  P E "<< potential_energy(w2,h);
      cout << "  E "<< kinetic_finite_diff(w2,h)+potential_energy(w2,h);
      cout << "  kin "<< q_expected_kinetic(w2,h);
      
      
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
	      double dipole_x=charge_nuclei*pow(soft_core+h.x1[i]*h.x1[i]+h.x2[j]*h.x2[j],-3./2.)*h.x1[i];
	      //	      double dipole_y=charge_nuclei*pow(soft_core+(h.x1[i])*(h.x1[i])+h.x2[j]*h.x2[j],-3./2.)*h.x2[j];
	      
	      double norm_temp=(w1.w[w1.index(k,j,i)][0]*w1.w[w1.index(k,j,i)][0]+w1.w[w1.index(k,j,i)][1]*w1.w[w1.index(k,j,i)][1]);
	      
	      dipole_x_axis+=w1.dx1*w1.dx2*w1.dx3*(dipole_x)*norm_temp;
	      
	    }

      fprintf(out7,"%e %e \n",abs(t),dipole_x_axis);
	  

      absorber(w1,0.1 ,0.1 ,0. ,0.1 ,0.1 ,0. ,1./6.);

    }//end time loop
  
}//End
