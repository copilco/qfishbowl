/*
 *  pulse.cpp
 *  
 *
 *  Created by Camilo Ruiz Mendez on 28/04/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "timegrid.h"
#include "timeobject.h"


class laser
{
	
public:
	
	double dt;
	double I0;
	double w;
	int cycles; /*Number of cycles in the cos2 pulse */
	double cep; /* Carrier envelope phase */

	double E0;
	double period;
	double lambda;
	int nmaxt;
	
	timegrid g;
	timeobject e,av;
	
	laser (double _dt, int _cycles, double _I0);
	laser(double _dt, double _I0, double _w, int _cycles, double _cep);

	
	void create_laser(double _dt, double _I0, double _w, int _cycles, double _cep )
	{
		dt=_dt;
		I0=_I0;
		w=_w;
		cycles=_cycles;
		cep=cep;
		
		E0=sqrt(I0/3.5e16);
		period=dospi/w;
		nmaxt=int(period*cycles/dt);
		
		g.set_grid(nmaxt,abs(dt));
		e.put_on_grid(g);
		av.put_on_grid(g);
		
		for(int kt=0;kt<e.n;kt++)
		{
			double t=e.t[kt];
			e.f[kt][0] =E0*cos(w*t+cep)*cos(w*t/2./cycles)*cos(w*t/2./cycles);
			e.f[kt][1]=0.;
			//Make a copy in the a vector which is later integrated.
			av.f[kt][0]=e.f[kt][0];
			av.f[kt][1]=e.f[kt][1];
		}
		av.integrate();
		
		double a=av.f[0][0];
		for(int kt=0;kt<e.n;kt++)
		{
			av.f[kt][0]-=a;
			av.f[kt][1]=0;			    
		}
		
		
	}//End create_laser;
		

	
};

laser::laser (double _dt, int _cycles, double _I0) {
	dt=_dt;
	I0=_I0;
	w=0.057;  //laser 800
	cycles=_cycles;
	cep=0.;
	
	E0=sqrt(I0/3.5e16);
	period=dospi/w;
	nmaxt=int(period*cycles/dt);
	
	g.set_grid(nmaxt,abs(dt));
	e.put_on_grid(g);
	av.put_on_grid(g);
	
	for(int kt=0;kt<e.n;kt++)
	{
		double t=e.t[kt];
		e.f[kt][0] =E0*cos(w*t+cep)*cos(w*t/2./cycles)*cos(w*t/2./cycles);
		e.f[kt][1]=0.;
		//Make a copy in the a vector which is later integrated.
		av.f[kt][0]=e.f[kt][0];
		av.f[kt][1]=e.f[kt][1];
	}
	av.integrate();
	
	double a=av.f[0][0];
	for(int kt=0;kt<e.n;kt++)
	{
		av.f[kt][0]-=a;
		av.f[kt][1]=0;			    
	}
	
}

laser::laser(double _dt, double _I0, double _w, int _cycles, double _cep) 
{
	dt=_dt;
	I0=_I0;
	w=_w;
	cycles=_cycles;
	cep=_cep;
	
	E0=sqrt(I0/3.5e16);
	period=dospi/w;
	nmaxt=int(period*cycles/dt);
	
	g.set_grid(nmaxt,abs(dt));
	e.put_on_grid(g);
	av.put_on_grid(g);
	
	for(int kt=0;kt<e.n;kt++)
	{
		double t=e.t[kt];
		e.f[kt][0] =E0*cos(w*t+cep)*cos(w*t/2./cycles)*cos(w*t/2./cycles);
		e.f[kt][1]=0.;
		//Make a copy in the a vector which is later integrated.
		av.f[kt][0]=e.f[kt][0];
		av.f[kt][1]=e.f[kt][1];
	}
	av.integrate();
	
	double a=av.f[0][0];
	for(int kt=0;kt<e.n;kt++)
	{
		av.f[kt][0]-=a;
		av.f[kt][1]=0;			    
	}
	
}//End create_laser;
	
	
	
