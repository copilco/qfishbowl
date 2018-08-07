#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "cartesian.h"


void placeWF(wavefunction wlarge, wavefunction wsmall)
{
  int kplacer=floor((wlarge.n3-wsmall.n3)/2);
  int jplacer=floor((wlarge.n2-wsmall.n2)/2);
  int iplacer=floor((wlarge.n1-wsmall.n1)/2);
  cout << "\n kplacer "<< kplacer
       << " jplacer "<< jplacer
       << " jplacer "<< iplacer;

for(int k=0;k<wlarge.n3*wlarge.n2*wlarge.n1;k++)
  { 
    wlarge.w[k][0]=0.;
    wlarge.w[k][1]=0.;	  
  }

  for(int k=0;k<wsmall.n3;k++)
    for(int j=0;j<wsmall.n2;j++)
      for(int i=0;i<wsmall.n1;i++)
	{ 	  
	  wlarge.w[wlarge.index(kplacer+k,jplacer+j,iplacer+i)][0]
	    =wsmall.w[wsmall.index(k,j,i)][0];
	  wlarge.w[wlarge.index(kplacer+k,jplacer+j,iplacer+i)][1]
	    =wsmall.w[wsmall.index(k,j,i)][1];
	}
}


void gaussian(wavefunction w)
{

  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{ 	  
	  w.w[w.index(k,j,i)][0]=   exp( -w.x1[i]*w.x1[i]-w.x2[j]*w.x2[j]-w.x3[k]*w.x3[k] )*sin(2.2*w.x2[j]);
	  //*cos(2.2*w.x2[j]);
	  w.w[w.index(k,j,i)][1]=0.;//  exp( -w.x1[i]*w.x1[i]-w.x2[j]*w.x2[j]-w.x3[k]*w.x3[k] )*sin(2.2*w.x2[j]);;
	}
}

void gaussianR0(wavefunction w,double x1,double x2,double x3)
{

  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{ 	  
	  w.w[w.index(k,j,i)][0]= exp( -(w.x1[i]-x1)*(w.x1[i]-x1)-(w.x2[j]-x2)*(w.x2[j]-x2)-(w.x3[k]-x3)*(w.x3[k]-x3) );//*cos(2.2*w.x2[j]);
	  w.w[w.index(k,j,i)][1]=0.;//  exp( -w.x1[i]*w.x1[i]-w.x2[j]*w.x2[j]-w.x3[k]*w.x3[k] )*sin(2.2*w.x2[j]);;
	}
}

void gaussianR1(wavefunction w, double kx1,double kx2,double kx3, double x01,double x02,double x03 )
{
  double rhox1=dospi/kx1;
  double rhox2=dospi/kx2;
  double rhox3=dospi/kx3;

  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
        {         
          w.w[w.index(k,j,i)][0]= exp( -(w.x1[i]-x01)*(w.x1[i]-x01)/rhox1/rhox1 - (w.x2[j]-x02)*(w.x2[j]-x02)/rhox2/rhox2
				       - (w.x3[k]-x03)*(w.x3[k]-x03)/rhox3/rhox3 )*cos(kx1*w.x1[i] + kx2*w.x2[j] + kx3*w.x3[k] );
          w.w[w.index(k,j,i)][1]= 0.;//exp( -(w.x1[i]-x01)*(w.x1[i]-x01)-(w.x2[j]-x02)*(w.x2[j]-x02)-(w.x3[k]-x03)*(w.x3[k]-x03) )*sin(4*w.x1[j]);
        }
}




void prop_kinetic(wavefunction &w, hamiltonian &h, complex dt)
{
  w.expected_kin=0.;
  w.qnorm=0.;
  double factor=1./w.n1/w.n2/w.n3;

  /**********************/
  fftw_execute(w.p3DF);
  /**********************/
  
  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{
	  
	  double kinetic=
	    w.q1[i]*w.q1[i]*h.a1+
	    w.q2[j]*w.q2[j]*h.a2+
	    w.q3[k]*w.q3[k]*h.a3;
	  
	  w.expected_kin+=w.dq1*w.dq2*w.dq3*kinetic*
	    (w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1])
	    *w.kfactor1*w.kfactor1*w.kfactor2*w.kfactor2*w.kfactor3*w.kfactor3;
	  
	  //The factor sould be: //(dx*dy*dz*dx*dy*dz/dospi/dospi/dospi);
	  
	    w.qnorm+=w.dq1*w.dq2*w.dq3*
	      (w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1])
	      *w.kfactor1*w.kfactor1*w.kfactor2*w.kfactor2*w.kfactor3*w.kfactor3;
	  
	  complex kinoperator=exp(-I*dt*kinetic );
	  
	  double aux0r=w.w[w.index(k,j,i)][0];
	  double aux0i=w.w[w.index(k,j,i)][1];
	  
	  w.w[w.index(k,j,i)][0]=( aux0r*real(kinoperator)+aux0i*imag(kinoperator) )*factor;
	  w.w[w.index(k,j,i)][1]=( aux0i*real(kinoperator)-aux0r*imag(kinoperator) )*factor;
	  
	  //w.w[w.index(k,j,i)][0]*=factor;
	  //w.w[w.index(k,j,i)][1]*=factor;
	  

	}
       
       /**********************/
       fftw_execute(w.p3DB);
       /**********************/


     
}

void prop_kinetic_laser_AP(wavefunction &w, hamiltonian &h, complex dt, double vect_pot1, double vect_pot2, double vect_pot3)
{

  w.expected_kin=0.;
  w.qnorm=0.;

  double factor=1./w.n1/w.n2/w.n3;
  /**********************/
  fftw_execute(w.p3DF);
  /**********************/
  
  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{
	  
	  double kinetic=
	    (w.q1[i]-h.b1*vect_pot1/lightC_au)*(w.q1[i]-h.b1*vect_pot1/lightC_au)*h.a1+
	    (w.q2[j]-h.b2*vect_pot2/lightC_au)*(w.q2[j]-h.b2*vect_pot2/lightC_au)*h.a2+
	    (w.q3[k]-h.b3*vect_pot3/lightC_au)*(w.q3[k]-h.b3*vect_pot3/lightC_au)*h.a3;
	  
	  w.expected_kin+=w.dq1*w.dq2*w.dq3*kinetic*
	    (w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1])
	    *w.kfactor1*w.kfactor1*w.kfactor2*w.kfactor2*w.kfactor3*w.kfactor3;
	  
	  //The factor sould be: //(dx*dy*dz*dx*dy*dz/dospi/dospi/dospi);
	  
	  w.qnorm+=w.dq1*w.dq2*w.dq3*
	    (w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1])
	    *w.kfactor1*w.kfactor1*w.kfactor2*w.kfactor2*w.kfactor3*w.kfactor3;


	  
	  complex kinoperator=exp(-I*dt*kinetic );
	  
	  double aux0r=w.w[w.index(k,j,i)][0];
	  double aux0i=w.w[w.index(k,j,i)][1];
	  
	  w.w[w.index(k,j,i)][0]=( aux0r*real(kinoperator)+aux0i*imag(kinoperator) )*factor;
	  w.w[w.index(k,j,i)][1]=( aux0i*real(kinoperator)-aux0r*imag(kinoperator) )*factor;
	  
	}
       
       /**********************/
       fftw_execute(w.p3DB);
       /**********************/

     
}

void prop_potential(wavefunction &w, hamiltonian h, complex dt )
{

  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{
	  complex potoperator=exp(-I*dt*h.v[h.index(k,j,i)]);
	  
	  double aux0r=w.w[w.index(k,j,i)][0];
	  double aux0i=w.w[w.index(k,j,i)][1];
	  
	  w.w[w.index(k,j,i)][0]=aux0r*real(potoperator)+aux0i*imag(potoperator);
	  w.w[w.index(k,j,i)][1]=aux0i*real(potoperator)-aux0r*imag(potoperator);
	}
  
}

void absorber(wavefunction w, double _frac_x1_left,double _frac_x2_left,double _frac_x3_left,double _frac_x1_right,double _frac_x2_right, double _frac_x3_right, double _exponent)
{
  
  
  double frac_x1_right=_frac_x1_right;
  double frac_x2_right=_frac_x2_right;
  double frac_x3_right=_frac_x3_right;
  
  double frac_x1_left=_frac_x1_left;
  double frac_x2_left=_frac_x2_left;
  double frac_x3_left=_frac_x3_left;
  
  double exponent=_exponent; //try use 1./6.

  double mask_start_x1_right=w.x1[int(w.n1*(1.-frac_x1_right))];
  double mask_start_x2_right=w.x2[int(w.n2*(1.-frac_x2_right))];
  double mask_start_x3_right=w.x3[int(w.n3*(1.-frac_x3_right))];
  
  double mask_start_x1_left=w.x1[int(w.n1*frac_x1_left)+1];
  double mask_start_x2_left=w.x2[int(w.n2*frac_x2_left)+1];
  double mask_start_x3_left=w.x3[int(w.n3*frac_x3_left)+1];
  
  double argument_x1_right;
  double argument_x2_right;
  double argument_x3_right;
  
  double argument_x1_left;
  double argument_x2_left;	
  double argument_x3_left;
  
  double mask_x1_right;
  double mask_x2_right;
  double mask_x3_right;
  
  double mask_x1_left;
  double mask_x2_left;
  double mask_x3_left;
  
  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{
	  argument_x1_right=(pi/2.)*(w.x1[i]-mask_start_x1_right)/(w.x1[w.n1-1]-mask_start_x1_right+1.e-20);
	  argument_x2_right=(pi/2.)*(w.x2[j]-mask_start_x2_right)/(w.x2[w.n2-1]-mask_start_x2_right+1.e-20);
	  argument_x3_right=(pi/2.)*(w.x3[k]-mask_start_x3_right)/(w.x3[w.n3-1]-mask_start_x3_right+1.e-20);
	  
	  argument_x1_left=(pi/2.)*(w.x1[i]-mask_start_x1_left)/(w.x1[0]-mask_start_x1_left+1.e-20);
	  argument_x2_left=(pi/2.)*(w.x2[j]-mask_start_x2_left)/(w.x2[0]-mask_start_x2_left+1.e-20);
	  argument_x3_left=(pi/2.)*(w.x3[k]-mask_start_x3_left)/(w.x3[0]-mask_start_x3_left+1.e-20);
	  
	  mask_x1_right=pow(fabs(cos(argument_x1_right)),exponent);
	  mask_x2_right=pow(fabs(cos(argument_x2_right)),exponent);
	  mask_x3_right=pow(fabs(cos(argument_x3_right)),exponent);
	  
	  mask_x1_left=pow(fabs(cos(argument_x1_left)),exponent);
	  mask_x2_left=pow(fabs(cos(argument_x2_left)),exponent);
	  mask_x3_left=pow(fabs(cos(argument_x3_left)),exponent);
	  
	  if (i< int(w.n1*frac_x1_left))
	    {		  
	      w.w[w.index(k,j,i)][0]*=mask_x1_left;	
	      w.w[w.index(k,j,i)][1]*=mask_x1_left;	
	    }
	  
	  if (i> int(w.n1*(1.-frac_x1_right)))
	    {		  
	      w.w[w.index(k,j,i)][0]*=mask_x1_right;	
	      w.w[w.index(k,j,i)][1]*=mask_x1_right;	
	    }
	  
	  if (j< int(w.n2*frac_x2_left))
	    {		  
	      w.w[w.index(k,j,i)][0]*=mask_x2_left;	
	      w.w[w.index(k,j,i)][1]*=mask_x2_left;	
	    }
	  
	  if (j> int(w.n2*(1.-frac_x2_right)))
	    {		  
	      w.w[w.index(k,j,i)][0]*=mask_x2_left;	
	      w.w[w.index(k,j,i)][1]*=mask_x2_left;	
	    }
	  
	  if (k< int(w.n3*frac_x3_left))
	    {		  
	      w.w[w.index(k,j,i)][0]*=mask_x3_left;	
	      w.w[w.index(k,j,i)][1]*=mask_x3_left;	
	    }
	  
	  if (k> int(w.n3*(1.-frac_x3_right)))
	    {		  
	      w.w[w.index(k,j,i)][0]*=mask_x3_left;	
	      w.w[w.index(k,j,i)][1]*=mask_x3_left;	
	    }
	}//End the loop on kji
  
}

double kinetic_finite_diff(wavefunction w, hamiltonian h)
{

  complex kinetic2=0.;
  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{
	  complex phi=complex(w.w[w.index(k,j,i)][0],w.w[w.index(k,j,i)][1]);
	  complex phi_ip=complex(0.,0.);
	  complex phi_im=complex(0.,0.);
	  complex phi_jp=complex(0.,0.); 
	  complex phi_jm=complex(0.,0.);	     
	  complex phi_kp=complex(0.,0.); 
	  complex phi_km=complex(0.,0.);
	  
	  /******************************************************************/
	  if (i-1>=0)
	    phi_ip=complex(w.w[w.index(k,j,i-1)][0],w.w[w.index(k,j,i-1)][1]);		
	  
	  if(i+1<w.n1)
	    phi_im=complex(w.w[w.index(k,j,i+1)][0],w.w[w.index(k,j,i+1)][1]);
	  
	  /******************************************************************/
	  if (j-1>=0)
	    phi_jp=complex(w.w[w.index(k,j-1,i)][0],w.w[w.index(k,j-1,i)][1]);
	  
	  if(j+1<w.n2)
	    phi_jm=complex(w.w[w.index(k,j+1,i)][0],w.w[w.index(k,j+1,i)][1]);
	  
	  /******************************************************************/
	  if (k-1>=0)
	    phi_kp=complex(w.w[w.index(k-1,j,i)][0],w.w[w.index(k-1,j,i)][1]);		
	  
	  if(k+1<w.n3)
	    phi_km=complex(w.w[w.index(k+1,j,i)][0],w.w[w.index(k+1,j,i)][1]);
	  
	  /******************************************************************/
	  
	  complex ksquare=
	    -h.a1*(phi_ip-2.*phi+phi_im)/w.dx1/w.dx1
	    -h.a2*(phi_jp-2.*phi+phi_jm)/w.dx2/w.dx2
	    -h.a3*(phi_kp-2.*phi+phi_km)/w.dx3/w.dx3;
	  
	  kinetic2+=w.dx1*w.dx2*w.dx3*(conj(phi)*ksquare);
	 
	}
  return real(kinetic2);
}

double potential_energy(wavefunction w, hamiltonian h)
{
  
  double potE=0.;

  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{
	  potE+=w.dx1*w.dx2*w.dx3*h.v[h.index(k,j,i)]*
	    (w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1]);
	}
  return potE;

}


void momentum_distribution(wavefunction &w , wavefunction &momentum_w , double x0, double y0, double z0)
{

  double r;
  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{
	  if(w.n3==1 && w.dx3==1.)
	    r=sqrt((w.x1[i]-x0)*(w.x1[i]-x0)+(w.x2[j]-y0)*(w.x2[j]-y0));
	  else
	    r=sqrt((w.x1[i]-x0)*(w.x1[i]-x0)+(w.x2[j]-y0)*(w.x2[j]-y0)+(w.x3[k]-z0)*(w.x3[k]-z0));
	  
	  double r0=10.;//5.;
	  double r1=20.;//15.;
	  
	  //double r0=20.;
	  //double r1=30.;
	  
	  if(r<=r0)
	    {
	      momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
	      momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
	    }
	  if(r>r0 && r<=r1 )
	    {
	      momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/16.));
	      momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/16.));;
	    }
	  if(r>r1)
	    {
	      momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
	      momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
	    }
	
      }

    fftw_execute(momentum_w.p3DF);

}


void Mask2D( wavefunction &w1 , double _r0, double _r1, double _sigma0 )
{
	
	double r0=_r0;//10.;//5.;
	double r1=_r1;//20.;//15.;
	
	double r;
	for(int k=0;k<w1.n3;k++)
		for(int j=0;j<w1.n2;j++)
			for(int i=0;i<w1.n1;i++)
			{
				r=sqrt( w1.x1[i]*w1.x1[i] + w1.x2[j]*w1.x2[j]  );
				
				if(r<=r0)
				{
					w1.w[w1.index(k,j,i)][0]*=0.;
					w1.w[w1.index(k,j,i)][1]*=0.;
				}
				if(r>r0 && r<=r1 )
				{
					w1.w[w1.index(k,j,i)][0]*=(1.-exp(-(r-r0)*(r-r0)/16.));
					w1.w[w1.index(k,j,i)][1]*=(1.-exp(-(r-r0)*(r-r0)/16.));;
				}
				if(r>r1)
				{
					w1.w[w1.index(k,j,i)][0]*=1.;
					w1.w[w1.index(k,j,i)][1]*=1.;
				}
								
			}
	
}


void AntiMask2D( wavefunction &w1 , double _r0, double _r1, double _sigma0 )
{
	
	double r0=_r0;//10.;//5.;
	double r1=_r1;//20.;//15.;
	
	double r;
	for(int k=0;k<w1.n3;k++)
		for(int j=0;j<w1.n2;j++)
			for(int i=0;i<w1.n1;i++)
			{
				r=sqrt( w1.x1[i]*w1.x1[i] + w1.x2[j]*w1.x2[j]  );
				
				/*
				if(r<=r0)
				{
					w1.w[w1.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
					w1.w[w1.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
					
				}*/
				
				if(r>r0 && r<=r1 )
				{
					w1.w[w1.index(k,j,i)][0]*=(exp(-(r-r0)*(r-r0)/_sigma0));
					w1.w[w1.index(k,j,i)][1]*=(exp(-(r-r0)*(r-r0)/_sigma0));;
				}
				if(r>r1)
				{
					w1.w[w1.index(k,j,i)][0]=0.;
					w1.w[w1.index(k,j,i)][1]=0.;
				}
				
			}
		
}

void momentum_distribution1D(wavefunction &w , wavefunction &momentum_w , double x0 )
{
	
	double r;
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
					r=sqrt( (w.x1[i]-x0)*(w.x1[i]-x0) );
				
				double r0=10.;//5.;
				double r1=30.;//15.;
				
				//double r0=20.;
				//double r1=30.;
				
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/50.));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/50.));;
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
				}
				
			}
	
	//  fftw_execute(momentum_w.p3DF);
	
}


void mask_function1D(wavefunction &w , wavefunction &momentum_w , double _r0, double _r1, double _sigma )
{
	
	double r;
	for(int k=0;k<w.n3;k++)
		for(int j=0;j<w.n2;j++)
			for(int i=0;i<w.n1;i++)
			{
				r=sqrt( (w.x1[i])*(w.x1[i]) );
				
				double r0=_r0;//10.;//5.;
				double r1=_r1;//30.;//15.;
				
				//double r0=20.;
				//double r1=30.;
				
				if(r<=r0)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=0.;
					momentum_w.w[momentum_w.index(k,j,i)][1]=0.;
				}
				if(r>r0 && r<=r1 )
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0]*(1.-exp(-(r-r0)*(r-r0)/_sigma));
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1]*(1.-exp(-(r-r0)*(r-r0)/_sigma));;
				}
				if(r>r1)
				{
					momentum_w.w[momentum_w.index(k,j,i)][0]=w.w[w.index(k,j,i)][0];
					momentum_w.w[momentum_w.index(k,j,i)][1]=w.w[w.index(k,j,i)][1];
				}
				
			}
	
	//  fftw_execute(momentum_w.p3DF);
	
}



double q_expected_kinetic(wavefunction &w, hamiltonian &h)
{
  w.expected_kin=0.;
  w.qnorm=0.;
  double factor=1./w.n1/w.n2/w.n3;

  /**********************/
  fftw_execute(w.p3DF);
  /**********************/
  
  for(int k=0;k<w.n3;k++)
    for(int j=0;j<w.n2;j++)
      for(int i=0;i<w.n1;i++)
	{
	  
	  double kinetic=
	    w.q1[i]*w.q1[i]*h.a1+
	    w.q2[j]*w.q2[j]*h.a2+
	    w.q3[k]*w.q3[k]*h.a3;
	  
	  w.expected_kin+=w.dq1*w.dq2*w.dq3*kinetic*
	    (w.w[w.index(k,j,i)][0]*w.w[w.index(k,j,i)][0]+w.w[w.index(k,j,i)][1]*w.w[w.index(k,j,i)][1])
	    *w.kfactor1*w.kfactor1*w.kfactor2*w.kfactor2*w.kfactor3*w.kfactor3;
	  	  
	  w.w[w.index(k,j,i)][0]*=factor;
	  w.w[w.index(k,j,i)][1]*=factor;
	  
	}
  
  
  /**********************/
  fftw_execute(w.p3DB);
  /**********************/

  return w.expected_kin;
}

complex projection(wavefunction &w1, wavefunction &w2)
{

  complex projection =complex(0.,0.);
  for(int k=0;k<w1.n3;k++)
    for(int j=0;j<w1.n2;j++)
      for(int i=0;i<w1.n1;i++)
	{
	  complex wave1=complex(w1.w[w1.index(k,j,i)][0],w1.w[w1.index(k,j,i)][1]);
	  complex wave2=complex(w2.w[w2.index(k,j,i)][0],w2.w[w1.index(k,j,i)][1]);
	  
	  projection+=w1.dx1*w1.dx2*w1.dx3*conj(wave2)*wave1;
	}

  return projection;

}


void project_out(wavefunction &w1, wavefunction &w2)
{

  complex p1=projection(w1,w2);

  for(int k=0;k<w1.n3;k++)
    for(int j=0;j<w1.n2;j++)
      for(int i=0;i<w1.n1;i++)
	{
	  complex wave1=complex(w1.w[w1.index(k,j,i)][0],w1.w[w1.index(k,j,i)][1]);
	  complex wave2=complex(w2.w[w2.index(k,j,i)][0],w2.w[w1.index(k,j,i)][1]);
	  
	  complex w=wave1-p1*wave2;
	  w1.w[w1.index(k,j,i)][0]=real(w);
	  w1.w[w1.index(k,j,i)][1]=imag(w);

	}

}

void snapshot(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3)
{

  for(int k=0;k<w1.n3/skiper3;k++)
    for(int j=0;j<w1.n2/skiper2;j++)
      for(int i=0;i<w1.n1/skiper1;i++)
	{
	  
	  double norm= w1.dx1*w1.dx2*w1.dx3*(
					     w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*
					     w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]+
					     w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*
					     w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]);
	  fprintf(file,"%e\n",norm);
	}



}

void Qsnapshot(FILE *file,wavefunction &w1, int skiper1, int skiper2, int skiper3)
{

	for(int k=0;k<w1.n3/skiper3;k++)
		for(int j=0;j<w1.n2/skiper2;j++)
			for(int i=0;i<w1.n1/skiper1;i++)
			{
				double qnorm=w1.dq1*w1.dq2*w1.dq3*
				(w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]*
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][0]+
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1]*
				 w1.w[w1.index(k*skiper3,j*skiper2,i*skiper1)][1])
				*w1.kfactor1*w1.kfactor1*w1.kfactor2*w1.kfactor2*w1.kfactor3*w1.kfactor3;
				fprintf(file,"%e\n",qnorm);
			}
}

/*
void prop( wavefunction &w)
{
  //cout << "  obs_kin "<< w.qnorm;
  w.expected_kin=6.;
  w.qnorm+=77.;
  //cout << "  obs_kin "<< w.qnorm;
}
*/
	


  
