/*
 *  pulse.cpp laser pulses train differents characteristics by pulse
 *  
 *
 *  Created by Alexis Chacón and Rafael Moran on 24/05/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef LASER03_H
#define LASER03_H
#include "timegrid.h"
#include "timeobject.h"
#include "MatrixC.h"
#include<vector>
using namespace std;

class laser03 
{

public:

   int nmaxt;
   int npulses;


   double dt;
   double t00;
   double alphaomega;


   vector<double> t;
   vector<double> I0;


   vector<double> w0;
   vector<double> cycles0;


   vector<double> cep0;
   vector<double> delay0;


   timegrid g;
   timeobject e,av;	


   /**************************************/
   /*  Creator Objet Laser and function
   /*************************************/
   laser03(int _npulses, int _nmaxt);
   void laser_pulses( double _dt, double _t0 );   


   /*********************************/
   /*  Mayor and minus functions
   /********************************/
   double mayor(vector<double>& v); 
   double minus(vector<double>& v);   

};

//Initialize variable
laser03::laser03(int _npulses, int _nmaxt )
{
   npulses=_npulses;
   nmaxt=_nmaxt;

   MatrixC fieldm(npulses, nmaxt);

   I0.resize(npulses,0.);
   w0.resize(npulses,0.);

   cycles0.resize(npulses,0.); 
   cep0.resize(npulses,0.);


   if (npulses ==1 )
      delay0.resize(npulses,0.);
   else 
      delay0.resize(npulses-1,0.);
}//End initialize variable  



   /*******************************/
   /*  START LASER FUNCTION
   /******************************/   
void laser03::laser_pulses( double _dt, double _t0 )
{
   MatrixC fieldm(npulses, nmaxt);

   dt=_dt;
   t00=_t0;

   //Start one pulse
   if (npulses==1)
   {
      double E00; 
      double w00;

      double cycles00; 
      double cep00;

      double period00; 
      double delay00; 

      double t00;
      double tm;
      E00=sqrt(I0[0]/3.5e16);

      w00=w0[0];
      cycles00=cycles0[0];

      cep00=cep0[0];
      period00=dospi/w00;

      delay00=delay0[0];	    
      tm=t00+cycles00*period00;

      alphaomega=period00;
      t.resize(nmaxt,0.);
	    	
      //axis time
      for(int k=0;k<nmaxt;k++)
 	 t[k]=t00-period00+k*dt;
		
     // objets of definition //
     g.set_grid(nmaxt,abs(dt));   
     e.put_on_grid(g);
     av.put_on_grid(g);

     //start loop creater pulse train
     for(int ki=0;ki<npulses;ki++)
     {	
	for(int kt=0;kt<nmaxt;kt++)
	{
 	   if (t[kt] < t00 || t[kt] > tm)        		
	        fieldm.Set(ki,kt,0.0,0.0);        		
	   else 
                fieldm.Set(ki,kt,
		       E00*sin(w00*(t[kt]-t00)/2/cycles00)*sin(w00*(t[kt]-t00)/2/cycles00)*sin(w00*(t[kt]-t00)+cep00),0.0);			
	}
     }//end loop creater pulse train
	

    //start loop sum pulses
    for (int kt=0;kt<nmaxt;kt++)
    {
	e.f[kt][0]=0;
	e.f[kt][1]=0;
		
    	for (int ki=0;ki<npulses;ki++)
	{
	   e.f[kt][0]+=fieldm.Get(ki,kt).r;
	   e.f[kt][1]+=fieldm.Get(ki,kt).i;  
	}	
	av.f[kt][0]=e.f[kt][0];
		   av.f[kt][1]=e.f[kt][1];		
    }//end loop sum pulses
	
    av.integrate();
    double a=av.f[0][0];
		
    for(int kt=0;kt<nmaxt;kt++)
    {
        av.f[kt][0]-=a; 
	av.f[kt][1]=0;    
    }	
   }//End one pulse

   //Start two or more pulses
   else
   {
      vector<double> E00; 
      vector<double> w00; 
      vector<double> cycles00;
      vector<double> cep00;
      vector<double> period00;
      vector<double> delay00;
      vector<double> t00;
      vector<double> tm;

      E00.resize(npulses,0.);
      w00.resize(npulses,0.);

      cycles00.resize(npulses,0.); 
      cep00.resize(npulses,0.);

      period00.resize(npulses,0.);
      delay00.resize(npulses-1,0.); 

      t00.resize(npulses,0.); 
      tm.resize(npulses,0.);


      for (int i=0;i<npulses;i++)
      {
	 E00[i]=sqrt(I0[i]/3.5e16); 
	 w00[i]=w0[i];

	 cycles00[i]= cycles0[i];
	 cep00[i]= cep0[i];

	 period00[i]= dospi/w00[i];
	
	 if (i<(npulses-1)) delay00[i]=delay0[i];
      }


      for (int k=0;k<npulses;k++)
      {		
	 if (k==0)
	 {
	    t00[k]=_t0;
	    tm[k]=t00[k]+cycles00[k]*period00[k];
	 }
	 else
         {
	    t00[k]=tm[k-1]+delay00[k-1]-(cycles00[k-1]*period00[k-1]+cycles00[k]*period00[k])/2;
	    tm[k]=t00[k]+cycles00[k]*period00[k];
	 } 
      }
		

      //Parametro de extremos
      alphaomega=period00[0];	
      for (int k=1;k<npulses;k++)	
	  alphaomega=alphaomega+period00[k];
      alphaomega=alphaomega/npulses;	

      double may=mayor(tm);
      double min=minus(t00);
      t.resize(nmaxt,0.);	

       //axis time
      for (int k=0;k<nmaxt;k++)
	  t[k]=min-alphaomega+k*dt;
	
      // objets of definition
      g.set_grid(nmaxt,abs(dt));
      e.put_on_grid(g);
      av.put_on_grid(g);
	
      //start loop pulse train differents
      for(int ki=0;ki<npulses;ki++)
      {	
	 for(int kt=0;kt<nmaxt;kt++)
	 {
	    if (t[kt] < t00[ki] || t[kt] > tm[ki])
                fieldm.Set(ki,kt,0.0,0.0);
        	    
	    else
	        fieldm.Set(ki,kt,
			E00[ki]*sin(w00[ki]*(t[kt]-t00[ki])/2/cycles00[ki])*sin(w00[ki]*(t[kt]-t00[ki])/2/cycles00[ki])*sin(w00[ki]*(t[kt]-t00[ki])+cep00[ki]),0.0);
	 }
      }//end loop pulse train differents 


      //start loop sum pulses
      for (int kt=0;kt<nmaxt;kt++)
      {
    	 e.f[kt][0]=0;
	 e.f[kt][1]=0;
		
    	 for (int ki=0;ki<npulses;ki++)
	 {
	    e.f[kt][0]+=fieldm.Get(ki,kt).r;
	    e.f[kt][1]+=fieldm.Get(ki,kt).i;
         }	
	 av.f[kt][0]=e.f[kt][0];
	 av.f[kt][1]=e.f[kt][1];		
      }
	
      av.integrate();
      double a=av.f[0][0];
		
      for(int kt=0;kt<nmaxt;kt++)
      {
           av.f[kt][0]-=a; 
	   av.f[kt][1]=0;			    
      }

   }//End two or more pulses
}//End laser_pulses


//Start function mayor
double laser03::mayor(vector<double>& v)
{
   int n=v.size();
   double may = v[0];
	
   for (int k=1;k<n;k++)   
       if (v[k]>may) may=v[k];
   
	return may;
}//End mayor


//Start function minus
double laser03::minus(vector<double>& v)
{
   int n=v.size();
   double min=v[0];
	

   for (int k=1;k<n;k++)
	if (v[k]<min) min=v[k];	

	return min;
}//End minus

#endif
//END
