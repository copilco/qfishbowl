/*
 *  laser.cpp laser pulses train differents characteristics by pulse
 *
 *  Created by Alexis Chacón and Camilo Ruiz
 *
 */

#ifndef LASER_H
#define LASER_H
#include <stdlib.h>
#include <string>
#include "timegrid.h"
#include "timeobject.h"
#include "grid.h"
#include "constant.h"
#include <vector>
using namespace std;

class laser 
{

public:
	
	string check;
	int Npulses;				//Total pulses number
	int Nmaxt;					//Total iteration number

	double t01;					//Initial time to first pulse
	double dt;					//Time step
	double blaser;				//Time before of the pulses
	double alaser;				//Time after of the pulses

	double major0;				//Major value
	double minus0;				//Minus value
	double max_twidth;
	
	
	vector<double> clock0;      //Initial time per pulse
	vector<double> clock1;      //The time until half per pulse
	vector<double> clock2;		//The time until final per pulse
	vector<double> twidth;		//Time bandwidth per pulse
	
	vector<int> kclock0;		//k-th value to the initial time per pulse 
	vector<int> kclock1;		//k-th value to the time until half per pulse
	vector<int> kclock2;		//k-th value to the time until final per pulse

	vector<int> Azeros;			//k-th zeros on the vector potential
	vector<int> Ezeros;			
	
	
	vector<double> I0;			//Intensity per pulse
	vector<double> e;			//Ellipticity per pulse
	vector<double> E0x;			//Electric Field Component in the x-direction
	vector<double> E0y;			//Electric Field Component in the y-direction
	
	
	vector<double> w0;			//Frequency per pulse
	vector<double> period0;		//Period per pulse
	vector<double> cycles0;		//Cycles number per pulse
	vector<double> chirp; 
	
	
	vector<double> cep0;		//Carrier envelope phase per pulse
	vector<double> phi_rel;		//Relative phase between Ex and Ey components of the electric field of the laser pulse
	vector<double> delay0;		//Delay between consecutive pulses
	
	
	vector<string> envelope;	//Envelope name used
	
	timegrid g;				//Time grid
	timeobject efield;		//Electric field for total pulses
	timeobject avector;		//Vector potential for total pulses
	
	timeobject *ef;			//Electric field per pulse
	timeobject *env;			//Envelope per pulse
	timeobject *av;			//Vector potential per pulse
	
	timeobject *av_int;		//Integral of the vector potential per pulse
	timeobject *avsq_int;	//Square integral of the vector potential per pulse	
	
	double *phase1;
	double *phase2; 
	double *induced_chirped_laser1;
	double *induced_chirped_laser2;	
	
	/*==========================*/
	/*      MAIN FUNCTIONS
	/*==========================*/
	laser(int _Npulses);      							                                // Creator Object
	//~laser();																			//Destructor
	void laser_pulses(double _dt, double _t01, double _blaser, double _alaser);	        // LASER PULSES TRAIN
	
	/*=========================*/
	/*   SECUNDARY FUNCTIONS
	/*=========================*/
	inline void Initialize_Amplitude_Period();			// Initialize max amplitude and period
	inline void Set_StartTime_EndTime();		        // "Start" and "end" time per pulse
	void Laser_Grid();									// Object time axis and field_set
	void Evaluation_Laser_Pulse();						// Evaluation pulses train
	void LSum_Pulses();									// Pulses sum
	void Set_Vector_Potential();						// Vector potential per pulse
	void Vector_Potential();							// Total vector potential
	void put_on_envelope();								// Generator of Envelope of the pulses
	void Set_Av_Integral();
	inline void Set_index_pulses();	
	void checking_starting_axis_time();
	
	double Lqmajor(vector<double>& v);					//Finding the maximum value of a vector
	double Lqminus(vector<double>& v);					//Finding the minimum value ot a vector
	
	
	void induced_chirp_laser(grid &space_grid, double Ip, double sw);
	void laser_coulomb_coupling(grid &space_grid, double ds, double xF);
	
	void laser_coulomb_couplingI(grid &space_grid, double ds, double xF);	
	void laser_coulomb_couplingII(grid &space_grid, double ds, double xF);
	
	void laser_coulomb_coupling_fxI(  grid &space_grid, double ds, double x01, double xF);  
	void laser_coulomb_coupling_fxII( grid &space_grid, double ds, double x02, double xF); 
	
	void laser_coulomb_coupling_x0(double *x0, grid &space_grid, double ds, double xF);
	void laser_coulomb_coupling_full(double xF, double ds, double soft_core, double Zatomic, grid &space_grid);	
	
	void laser_coulomb_coupling_x0I(grid &space_grid, double ds, double x01, double xF);
	void laser_coulomb_coupling_x0II(grid &space_grid, double ds, double x01, double xF);	
	
	void zeros_electric_f();
	void zeros_vector_pot();
	int vector_potential_match(double tI,double tII);
	
	void pulses_display();
	void saving_pulses(FILE *timeaxis,FILE *file, FILE *integrals);

	
	void save_ezeros(FILE *file);
	void save_azeros(FILE *file);
	
};


//==============================================================================//
			/*=== MAIN FUNCTIONS ===*/

/*===========================================================          		
		/*=== OBJECT'S CONSTRUCTOR LASER  ===*/
laser::laser(int _Npulses)
  {

	  
	  Npulses=_Npulses;
	  clock0.resize(Npulses,0.0);
	  clock1.resize(Npulses,0.0);	  
	  clock2.resize(Npulses,0.0);
	  
	  kclock0.resize(Npulses,0.0);
	  kclock1.resize(Npulses,0.0);
	  kclock2.resize(Npulses,0.0);
	  
	  twidth.resize(Npulses,0.0);
	  
	  I0.resize(Npulses,0.0);
	  e.resize(Npulses,0.0);	  
	  
	  envelope.resize(Npulses);
	  
	  for (int kpulse=0; kpulse<Npulses; kpulse++) 
		 envelope[kpulse] ="gauss";
	  
	  check="DEFAULT";
	  
	  
	  E0x.resize(Npulses,0.0);
	  E0y.resize(Npulses,0.0);	  
     
	  //Azeros.resize(Npulses);			
	  //Ezeros;	  
	  
	  w0.resize(Npulses,0.0);
	  period0.resize(Npulses,0.0);
     
	  
	  chirp.resize(Npulses,0.);
	  cycles0.resize(Npulses,0.0); 
	  cep0.resize(Npulses,0.0);	  
	  phi_rel.resize(Npulses,0.0);
	  
	  	  
	  ef  = new timeobject[Npulses];				// reserve memory
	  env = new timeobject[Npulses];				// reserve memory
	  av  = new timeobject[Npulses];	  
	  
	  
	  av_int    = new timeobject[Npulses];		// reserve memory
	  avsq_int  = new timeobject[Npulses];		// reserve memory	  	  
	  
	  if (Npulses > 1 )
		  delay0.resize(Npulses-1,0.0);	  		  
  }//End initialize variable  



/*laser::~laser(){
	delete avsq_int;
	delete av_int;
	delete ef;
	delete env;
	delete av;
}*/


/*===========================================================          		
 /*===  LASER PULSES FUNCTION  ===*/
void laser::laser_pulses(double _dt, double _t01, double _blaser, double _alaser)
{
	//return;
   	/*======= Pulse's Parameter ======*/
   	t01		= _t01;										//	Start time first pulse
   	dt		= abs(_dt);			                        //	Time step 	
   	
	blaser	= abs(_blaser);								//	Time before the laser or the pulse train
   	alaser	= abs(_alaser);								//	Time after the laser or the pulse train 	
	
	
   	//================================//
	Initialize_Amplitude_Period();						//	Initialize max amplitude and period        
   	Set_StartTime_EndTime();							//	"Start" and "end" time by each pulse
	
	
   	major0	= Lqmajor(clock2);	     	                //	Maximum time of all pulse
   	minus0	= Lqminus(clock0);							//	Minimun time of all pulse 
	
		
	//CHECKING CONTROL OF AXIS 
	checking_starting_axis_time();	
	Set_index_pulses();		
	
	
	
   	Laser_Grid();										//	Objet time axis      
   	Evaluation_Laser_Pulse();							//	Evaluation pulses
	
   	
	
	LSum_Pulses();										//	Pulses sum (build pulses train)
	Set_Vector_Potential();								//	Vector potential by each pulse
   	Vector_Potential();									//	Vector potential to pulses train 
	

	

}//End laser_pulses

/*============ END MAIN FUNCTIONS ===============*/



/*=================================================================*/
		/*========= SECUNDARY FUNCTIONS ===========*/

//== Function initialize max amplitude and period ==// 
inline void laser::Initialize_Amplitude_Period()
{
	
	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{
		double factor      =  1.0/(1.0+e[kpulse]*e[kpulse]);           //The ellipticity is e = E0x/E0y, Other condition
		E0x[kpulse]        =  e[kpulse]*sqrt(I0[kpulse]*factor/3.5e16);
		E0y[kpulse]        =  sqrt(I0[kpulse]*factor/3.5e16); 		
		period0[kpulse]    =  dospi/w0[kpulse];
		twidth[kpulse]	   =  cycles0[kpulse]*period0[kpulse];
	}
}


//== Function "start" and "end" time by each pulse ==//
inline void laser::Set_StartTime_EndTime()
{
	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{		
		if(kpulse==0)
		{
			clock0[kpulse]   = t01;
			clock1[kpulse]   = clock0[kpulse] + twidth[kpulse]/2.0;
			clock2[kpulse]   = clock0[kpulse] + twidth[kpulse];
		}
		else
		{

			clock0[kpulse]  = clock2[kpulse-1] + delay0[kpulse-1]
			                    - (twidth[kpulse-1] + twidth[kpulse])/2.0;
			
			clock1[kpulse]  = clock0[kpulse] + twidth[kpulse]/2.0;
			
			clock2[kpulse]  = clock0[kpulse] + twidth[kpulse];
		} 
	}//End loop
}//End function


inline void laser::Set_index_pulses()
{
	double min_aux0 = 0;

	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{		
		kclock0[kpulse] = floor( abs((clock0[kpulse]) - (minus0-blaser))/dt)	+ 1;
		kclock1[kpulse] = floor( abs((clock1[kpulse]) - (minus0-blaser))/dt )	+ 1;
		kclock2[kpulse] = floor( abs((clock2[kpulse]) - (minus0-blaser))/dt )	+ 1;			
	}//End loop

}

void laser::checking_starting_axis_time(){
	
	int s=0;
	if(check=="ControlGrid")
		s=1;

	max_twidth			=  twidth[0];//Lqmajor(twidth); 
	double start_time	=  t01 - blaser;			
	double final_time	=  t01+max_twidth+alaser;
	
	
	switch (s) {
		case 1:
			
			for (int kpulse=1; kpulse<Npulses; kpulse++) {
				
				if (start_time>minus0)
				{
					cout << "\n\n//***** ERROR FROM LASER PULSES: *******\n Please, check offset time or blaser parameter *****//\n";
					exit (1);
				}
				if(final_time<major0)			
				{
					cout << "\n\n//***** ERROR FROM LASER PULSES:  *******\n Please, check time after pulses or alaser parameter *****//\n";
					exit (1);
				}
			};			
	};
	
	minus0				= t01;
	major0				= final_time;

}


//== Function object time axis and field set  ==//
void laser::Laser_Grid()   
{	
	double final_time	= major0+alaser;
	double initial_time = minus0-blaser;
	Nmaxt				= floor((final_time - initial_time)/dt)+1;
	
	if (Nmaxt%2!=0) 
		Nmaxt=Nmaxt+1;	

	
	g.set_grid(Nmaxt, dt, initial_time);      
	efield.put_on_grid(g);
	avector.put_on_grid(g);
	
	for (int kfield=0; kfield <Npulses;kfield++){
		env[kfield].put_on_grid(g);
		ef[kfield].put_on_grid(g);
		av[kfield].put_on_grid(g);
	}
	cout << "\n/********************************/\n";
	cout << "        About laser.h ";	
	cout << "\nNtime= "<<Nmaxt;
	cout << "      Memory= "<<((2+3*Npulses)*24+8)*Nmaxt*1e-6 << "  Mb"<<endl;
	cout << "/********************************/\n";		
}


//== Funtion that evaluate of the pulses train ==//
void laser::Evaluation_Laser_Pulse()
{	
	put_on_envelope();
	
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{
				
		double a			= 64.*log10(2.)/twidth[kpulse]/twidth[kpulse];	
		double chirp_width  = 4.*a*chirp[kpulse];	
		
		
		for(int ktime=0;ktime<Nmaxt;ktime++)
		{
			
			double arg1=  w0[kpulse]*(g.t[ktime] - clock1[kpulse] )
					      + a*chirp_width*pow( g.t[ktime] - clock1[kpulse] ,2.)/(1. + chirp_width*chirp_width ) 
						  + cep0[kpulse]; 
			
			double arg2=  w0[kpulse]*(g.t[ktime] - clock1[kpulse] )
						  + a*chirp_width*pow( g.t[ktime] - clock1[kpulse] ,2.)/(1. + chirp_width*chirp_width )
						  + cep0[kpulse]+phi_rel[kpulse]; 
			
				
			ef[kpulse].f[ktime][0] = env[kpulse].f[ktime][0]*sin(arg1);			   
			ef[kpulse].f[ktime][1] = env[kpulse].f[ktime][1]*sin(arg2);
			
			
			av[kpulse].f[ktime][0]=ef[kpulse].f[ktime][0];
			av[kpulse].f[ktime][1]=ef[kpulse].f[ktime][1];		   
	   }
	}
}


//Place Laser Pulse Envelope 
void laser::put_on_envelope(){
	
	//string _name_env=envelope;
	vector<string> envelopes (5);
	
	envelopes[0]="rect";
	envelopes[1]="sin2";		
	envelopes[2]="rsin2";
	envelopes[3]="gauss";	
	envelopes[4]="konst";
	
	
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{	
	//Rectangle Envelope 
	if (envelope[kpulse]==envelopes[0]) {		
		

			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
				if (g.t[ktime] < clock0[kpulse] || g.t[ktime] > clock2[kpulse]){					
					env[kpulse].f[ktime][0] = 0.0;
					env[kpulse].f[ktime][1] = 0.0;
				} 
				else{
					env[kpulse].f[ktime][0] = E0x[kpulse];			   
					env[kpulse].f[ktime][1] = E0y[kpulse];			   
					
				}
				
			}
		}//End Rectangle Envelope
		
		

		//sin2 Envelope 
		if (envelope[kpulse]==envelopes[1]) 			
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				double arg0=w0[kpulse]*(g.t[ktime]-clock0[kpulse])/2./cycles0[kpulse];            //Envelope Argument
				env[kpulse].f[ktime][0] = E0x[kpulse]*sin(arg0)*sin(arg0);			   
				env[kpulse].f[ktime][1] = E0y[kpulse]*sin(arg0)*sin(arg0);			   
				
			}
		//End sin2 Envelope
		
		
		//Rectangle Multiplied sin2 Envelope (rsin2) 
		if (envelope[kpulse]==envelopes[2]) {		
			
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
				if (g.t[ktime] < clock0[kpulse] || g.t[ktime] > clock2[kpulse]){					
					env[kpulse].f[ktime][0] = 0.0;
					env[kpulse].f[ktime][1] = 0.0;
				} 
				else{
					double arg0=w0[kpulse]*(g.t[ktime]-clock0[kpulse])/2./cycles0[kpulse];            //Envelope Argument
					
					env[kpulse].f[ktime][0] = E0x[kpulse]*sin(arg0)*sin(arg0);			   
					env[kpulse].f[ktime][1] = E0y[kpulse]*sin(arg0)*sin(arg0); 
					
				}
				
			}			
		}//End Rectangle  Multiplied sin2 Envelop (rsin2)
		

		
		//Gauss Envelope 
		if (envelope[kpulse]==envelopes[3])
		{					
			double a			= 64.*log10(2.)/twidth[kpulse]/twidth[kpulse];
			double chirp_width	= 4.*a*chirp[kpulse];			
			
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
				env[kpulse].f[ktime][0] = E0x[kpulse]*
						exp( -a*( g.t[ktime] - clock1[kpulse] )*
							    ( g.t[ktime] - clock1[kpulse] )
								/(1. + chirp_width*chirp_width) ); 
				
				env[kpulse].f[ktime][1] = E0y[kpulse]*
						exp( -a*( g.t[ktime] - clock1[kpulse] )*
							    ( g.t[ktime] - clock1[kpulse] )
								/(1. + chirp_width*chirp_width) );
					
			}//End Gauss Envelope	
		
		}
		
		//Constant Envelope 
		if (envelope[kpulse]==envelopes[4]) {		
			
			for(int ktime=0;ktime<Nmaxt;ktime++)
			{
				
				env[kpulse].f[ktime][0] = E0x[kpulse];			   
				env[kpulse].f[ktime][1] = E0y[kpulse];			   
				
			}
						
		}//End Constant Envelope		
	}//End number of pulses
	
	
}



//== Function pulses sum (build pulses train) ==//
void laser::LSum_Pulses()
{
	//Start loop sum pulses
	for(int ktime=0;ktime<Nmaxt;ktime++)
	{
		efield.f[ktime][0]=0.0;
		efield.f[ktime][1]=0.0;		
		for(int kpulse=0;kpulse<Npulses;kpulse++)
		{
			efield.f[ktime][0]+= ef[kpulse].f[ktime][0];
			efield.f[ktime][1]+= ef[kpulse].f[ktime][1];

		}	
	}//End sum set pulses       
}//*/



//== Function vector potential by each pulse ==//
void laser::Set_Vector_Potential()
{
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{
		av[kpulse].integrateRK4();//av[kpulse].integrate();      //
		double aa=av[kpulse].f[0][0];
		double bb=av[kpulse].f[0][1];
		for(int ktime=0;ktime<Nmaxt;ktime++)
		{
			av[kpulse].f[ktime][0]-=aa; 
			av[kpulse].f[ktime][1]-=bb;
			
			av[kpulse].f[ktime][0]*=-1; 	//According to Landau & Lifshitz in The Classical Theory of Fields,
			av[kpulse].f[ktime][1]*=-1; 	//and Remetter Attosecond Wave Packet Interference, vector potential has negative sign.
			
			
		}  
	}
}





//== Function Integral of the vector potential by each pulse ==//
void laser::Set_Av_Integral()
{
	for(int kpulse=0;kpulse<Npulses;kpulse++)
	{		
		av_int[kpulse].put_on_grid(g);
		avsq_int[kpulse].put_on_grid(g);
		
		
		for(int ktime=0;ktime<Nmaxt;ktime++)
		{
			
			av_int[kpulse].f[ktime][0]=av[kpulse].f[ktime][0]; 
			av_int[kpulse].f[ktime][1]=av[kpulse].f[ktime][1];
			
			
			avsq_int[kpulse].f[ktime][0]=av[kpulse].f[ktime][0]*av[kpulse].f[ktime][0]; 
			avsq_int[kpulse].f[ktime][1]=av[kpulse].f[ktime][1]*av[kpulse].f[ktime][1];
			
		}
		
		av_int[kpulse].integrateRK4();		//av[kpulse].integrate(); 
		avsq_int[kpulse].integrateRK4();	//av[kpulse].integrate(); 				
		
	}
	
} //End integral of the vector potential





//== Function: vector potential to pulses train ==//
void laser::Vector_Potential()
{
	for(int ktime=0;ktime<Nmaxt;ktime++)
	{
		avector.f[ktime][0]=0.0;
		avector.f[ktime][1]=0.0;		
		for(int kpulse=0;kpulse<Npulses;kpulse++)
		{
			avector.f[ktime][0]+= av[kpulse].f[ktime][0];
			avector.f[ktime][1]+= av[kpulse].f[ktime][1];
		}	
	}//End sum set pulses       
}



//Start function major
double laser::Lqmajor(vector<double>& v)
{
	int n=v.size();
	double may = v[0];
	
	for (int k=1;k<n;k++)   
		if (v[k]>may) may=v[k];
	
	return may;
}//End major



//Start function minus
double laser::Lqminus(vector<double>& v)
{
	int n=v.size();
	double min=v[0];
	
	for (int k=1;k<n;k++)
		if (v[k]<min) min=v[k];	
	
	return min;
}//End minus






//Member function to calculate zeros on the electric field
void laser::zeros_electric_f()
{
	int Nzeros=2*cycles0[0]+2;
	Ezeros.resize(Nzeros, 0);
	
	
	int kont				= 1;		//kont
	int flag;							//flag
	
	
	int kinitial			= kclock0[0]+10;//floor(0.5*period0[0]/4./dt);
	int kfinal				= kclock2[0]-10;//-floor(0.5*period0[0]/4./dt);
	Ezeros[0]				= kclock0[0];
	Ezeros[Nzeros-1]		= kclock2[0];	

		
	double slope = (ef[0].f[kinitial+1][1]-ef[0].f[kinitial][1])/dt;
	
	//Finding good condition for flag
	if (slope<0.)
	{
		if(ef[0].f[kinitial+1][1]>0.)
			flag=0;
		else 
			flag=1;	
	}
	else{ 
		if(ef[0].f[kinitial+1][1]<0.)
			flag=1;
		else 
			flag=0;	
	}
	
	
	//Starting loop to find zeros
	for(int ktime=kinitial;ktime<kfinal;ktime++)
	{
		slope = (ef[0].f[ktime+1][1]-ef[0].f[ktime][1])/dt;		
		if(slope<0.){
			if ( (ef[0].f[ktime][1]<0.)&(flag==0) ) {
				Ezeros[kont] = ktime-1;
				kont+=		 1;
				flag		=1;
			}
		}
		else{ 
			if((ef[0].f[ktime][1]>0.)&(flag==1))
			{
				Ezeros[kont] = ktime-1;
				kont+=		   1;
				flag		 = 0;
			}
		}
		//Break out
		if (kont==Nzeros-1) 
			break;
		
	}//End loop to find zeros       	
}







//Member function to calculate zeros on laser potential vector
void laser::zeros_vector_pot()
{
	int Nzeros=2*cycles0[0]+1;	
	Azeros.resize(Nzeros, 0.);
		
	int kont	 = 1;
	int flag;
	
	int kinitial	 = kclock0[0]+floor(.5*period0[0]/4./dt);
	int kfinal		 = kclock2[0]-floor(.5*period0[0]/4./dt);
	
	Azeros[0]        = kclock0[0];
	Azeros[Nzeros-1] = kclock2[0];	
	
	//Slope 
	double slope = (av[0].f[kinitial+1][1]-av[0].f[kinitial][1])/dt;
	
	
	//Finding good condition for flag
	if (slope<0.)
	{
		if(av[0].f[kinitial][1]>0.)
			flag=0;
		else 
			flag=1;	
	}
	else{ 
		if(av[0].f[kinitial][1]<0.)
			flag=1;
		else 
			flag=0;	
	}
	
	
	//Starting loop to find zeros
	for(int ktime=kinitial;ktime<kfinal;ktime++)
	{
		slope = (av[0].f[ktime+1][1]-av[0].f[ktime][1])/dt;		
		if(slope<0.){
			if ( (av[0].f[ktime][1]<0.)&(flag==0) ) {
				Azeros[kont] = ktime-1;
				kont+=		 1;
				flag		=1;
			}
		}
		else{ 
			if((av[0].f[ktime][1]>0.)&(flag==1))
			{
				Azeros[kont] = ktime-1;
				kont+=		   1;
				flag		 = 0;
			}
		}
		
		//Break out
		if (kont==Nzeros-1) 
			break;		
		
	}//End loop to find zeros       	
}//End loop zeros in electric field





//Starting induced laser chirp Routine
void laser::induced_chirp_laser(grid &space_grid, double Ip, double sw)
{
	//Here I'm evaluating the chirped-IR-laser-induced 
	//For photoionization of on electron by the XUV attosecond pulses
	induced_chirped_laser1 =(double*)malloc(space_grid.n1*sizeof(double));	
	
	double Energy01;
	double alpha1; 
	double Energy02;
	double alpha2;	
	
	if (sw==0){
		Energy01 = w0[1]-Ip;
		alpha1   = 64.*log10(2.)/twidth[1]/twidth[1]; 
		
		//Starting chirped loop 
		for (int kmom=0; kmom<space_grid.n1; kmom++) {
			double betha_k = (space_grid.k1[kmom] + av[0].f[kclock1[1]][1] )*ef[0].f[kclock1[1]][1];
			double omega_k = pow(space_grid.k1[kmom] + av[0].f[kclock1[1]][1],2.)/2.-Energy01;
			double gamma_k = pow((4.*alpha1),2.) + pow((2.*betha_k),2.);
			
			induced_chirped_laser1[kmom]=2*betha_k*omega_k*omega_k/gamma_k;
			
		}//End loop evaluation of the chirped formula
	}
	else {
		induced_chirped_laser2=(double*)malloc(space_grid.n1*sizeof(double));	
		
		Energy01 = w0[1]-Ip;
		alpha1   = 64.*log10(2.)/twidth[1]/twidth[1]; 				
		
		Energy02 = w0[2] - Ip;
		alpha2   = 64.*log10(2.)/twidth[2]/twidth[2];
				
		
		//Starting chirped loop 
		for (int kmom=0; kmom<space_grid.n1; kmom++) {
			
			double betha_k = (space_grid.k1[kmom] + av[0].f[kclock1[1]][1] )*ef[0].f[kclock1[1]][1];
			double omega_k = pow(space_grid.k1[kmom] + av[0].f[kclock1[1]][1],2.)/2.-Energy01;
			double gamma_k = pow((4.*alpha1),2.) + pow((2.*betha_k),2.);
			
			induced_chirped_laser1[kmom]=-2*betha_k*omega_k*omega_k/gamma_k;			
			
			
			betha_k        = (space_grid.k1[kmom] + av[0].f[kclock1[2]][1] )*ef[0].f[kclock1[2]][1];
			omega_k        = pow(space_grid.k1[kmom] + av[0].f[kclock1[2]][1],2.)/2. - Energy02;
			gamma_k        = pow((4.*alpha2),2.) + pow((2.*betha_k),2.);
			
			induced_chirped_laser2[kmom]=2*betha_k*omega_k*omega_k/gamma_k;			
		}//End loop evaluation of the chirped formula
	}

}//End Routine to evaluate chirped formula the on time




//Laser-Coulomb Coupling IR
void laser::laser_coulomb_coupling(grid &space_grid, double ds, double xF)
{
	//Here I integrate Manfred's formula for
	//Laser-Coulomb Coupling between IR laser and coulomb potential well
	//The formula is integrated on time 
	
	
	phase1=(double*)malloc(space_grid.n1*sizeof(double));
	phase2=(double*)malloc(space_grid.n1*sizeof(double));
	
	
	double x01;
	double s;
	double aux_xF=xF;
	double integral_on_s;
	double integral_on_s2;
	double integral_on_time;
	double taux1;
	double taux2;	
	long double aux_index;
	
	int Ns;	
	int final_ktime;	

	//Set_Av_Integral();
	
	
	
	
	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {
		
		x01		= 1./(2.*(space_grid.k1[kmom] + av[0].f[kclock1[1]][1] ) );		
		
		
		if (space_grid.k1[kmom]<0){
			ds = -abs(ds);
			xF = -abs(aux_xF);
		}
		else{
			ds = +abs(ds);
			xF = +abs(aux_xF);
		}
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((xF-x01)/ds)) ;
		
		
		
		
		//Loop on s variable 
		integral_on_s			= 0.;
		integral_on_s2			= 0.;			
		for (int kspace=0; kspace<Ns; kspace++) {			
						
			//Axis of the variable s
			s					= x01 + kspace*ds;			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			taux1				= (g.t[kclock1[1]] + abs((s-x01)/(space_grid.k1[kmom])) - g.t[0]);
			aux_index			= floor(taux1/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[1]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[1]][1];	
			//End sentece on time (time integration)
			
			
			
			
			
			//double soft_core = 0.5;
			integral_on_s+=  - integral_on_time*ds/s/s
							 + integral_on_time*integral_on_time*pow(s,-3.)*ds+ 
							 - integral_on_time*integral_on_time*integral_on_time*pow(s,-4.)*ds;
			
			
			
			//Second part or second phase using the pulse
			x01 = 1./(2.*(space_grid.k1[kmom] + av[0].f[kclock1[2]][1] ) );		
			
			s					= x01 + kspace*ds;
			
			
			//Sentece on time
			integral_on_time	= 0.;
			final_ktime			= 0;
			
			taux2				= (g.t[kclock1[2]] + abs((s-x01)/(space_grid.k1[kmom])) - g.t[0]);
			aux_index			= floor(taux2/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time= av_int[0].f[Nmaxt-1][1]-av_int[0].f[kclock1[2]][1];					
			else 
				integral_on_time= av_int[0].f[final_ktime-1][1]-av_int[0].f[kclock1[2]][1];
			//End sentece on time
			
			
			
			
			integral_on_s2+=  - integral_on_time*ds/s/s + 
							  + integral_on_time*integral_on_time*pow(s,-3.)*ds +  
							  - integral_on_time*integral_on_time*integral_on_time*pow(s,-4.)*ds;
			
			
			
		}//End loop s variable
		phase1[kmom] = integral_on_s/abs(space_grid.k1[kmom]);
		phase2[kmom] = integral_on_s2/abs(space_grid.k1[kmom]);		
		
		
		
		if(kmom%(space_grid.n1/1)==0){
			cout << "\nIntegral1= " << integral_on_s;	
			cout << "     Integral2= " << integral_on_s2;			
		}		
		
	}//End momentum loop
	
}//End Manfred's Formula integration on time*/


//Laser-Coulomb Coupling IR
void laser::laser_coulomb_couplingI(grid &space_grid, double ds, double xF)
{
	//Here I integrate Manfred's formula for
	//Laser-Coulomb Coupling between IR laser and Coulomb potential well
	//The formula is integrated on time 
	
	
	phase1=(double*)malloc(space_grid.n1*sizeof(double));
	
	
	double x01;
	double s;
	double aux_xF=xF;
	double integral_on_s;
	double integral_on_time;
	double taux1;
	long double aux_index;
	
	int Ns;	
	int final_ktime;	
		
	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {
		
		x01		= 1./(2.*(space_grid.k1[kmom] + av[0].f[kclock1[1]][1] ) );		
		
		
		if (space_grid.k1[kmom]<0){
			ds = -abs(ds);
			xF = -abs(aux_xF);
		}
		else{
			ds = +abs(ds);
			xF = +abs(aux_xF);
		}
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((xF-x01)/ds)) ;
		
		
		
		//Loop on s variable 
		integral_on_s		= 0.0;
		for (int kspace=0; kspace<Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x01 + kspace*ds;			
			
			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			taux1				= (	g.t[kclock1[1]] + abs((s-x01)/space_grid.k1[kmom] ) - g.t[0]);
			aux_index			=	floor(taux1/dt);
			
			
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[1]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[1]][1];	
			//End sentece on time (time integration)
			
			
			double soft_core	= 0.0;
			integral_on_s+= integral_on_time*ds/s/s;//-pow( soft_core + pow( integral_on_time+s, 2. ),-1./2. )*ds;//integral_on_time*ds/s/s;//
			//+integral_on_time*integral_on_time*ds/s/s/s;
			
		}//End loop s variable
		phase1[kmom] = -integral_on_s/(space_grid.k1[kmom]);
			
	}//End momentum loop
	
}//End Manfred's Formula integration on time*/




void laser::laser_coulomb_couplingII(grid &space_grid, double ds, double xF)
{
	//Here I integrate Manfred's formula for
	//Laser-Coulomb Coupling between IR laser and coulomb potential well
	//The formula is integrated on time 
	
	
	phase2=(double*)malloc(space_grid.n1*sizeof(double));
	
	
	double x01;
	double s;
	double aux_xF	= xF;
	double integral_on_s;
	double integral_on_time;
	double taux1;
	long double aux_index;
	
	int Ns;	
	int final_ktime;	
	
	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {
		
		x01		= 1./(2.*(space_grid.k1[kmom] + av[0].f[kclock1[2]][1] ) );		
		
		
		if (space_grid.k1[kmom]<0){
			ds = -abs(ds);
			xF = -abs(aux_xF);
		}
		else{
			ds = +abs(ds);
			xF = +abs(aux_xF);
		}
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((xF-x01)/ds)); 
		
		
		
		
		//Loop on s variable 
		integral_on_s		= 0.;
		for (int kspace=0; kspace<Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x01 + kspace*ds;			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			taux1				= (g.t[kclock1[2]] + abs((s-x01)/(space_grid.k1[kmom])) - g.t[0]);
			aux_index			= floor(taux1/dt);
			
			
			if (int(aux_index)>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[2]][1];
			else
				integral_on_time = av_int[0].f[int(aux_index)-1][1] - av_int[0].f[kclock1[2]][1];	
			//End sentece on time (time integration)			
			
			
			
			
			double soft_core = 0.0;
			integral_on_s+= integral_on_time*ds/s/s;//-pow( soft_core +pow(integral_on_time+s,2.),-1./2.)*ds;
			//integral_on_time*ds/s/s;
			
		}//End loop s variable
		phase2[kmom] = -integral_on_s/(space_grid.k1[kmom]);
		
	}//End momentum loop
	
}//End Manfred's Formula integration on time*/





//Laser-Coulomb Coupling IR
void laser::laser_coulomb_coupling_x0I(grid &space_grid, double ds, double x01, double xF)
{
	//Here I integrate Manfred's formula for
	//Laser-Coulomb Coupling between IR laser and coulomb potential well
	//The formula is integrated on time 
	
	
	phase1=(double*)malloc(space_grid.n1*sizeof(double));
	
	
	double s;
	double aux_xF			= xF;
	
	
	double integral_on_s	;
	double integral_on_time	;
	
	
	double taux1;
	double x0_prime			= 0.;
	
	
	double soft_core		= 0.;
	long double aux_index;
	
	
	int Ns;	
	int final_ktime;	
	
	
	
	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {
		
		//Step space 
		if (space_grid.k1[kmom]<0){
			ds	= -abs(ds);
			xF	= -abs(aux_xF);
			x01	= -abs(x01);
		}
		else{
			ds	= abs(ds);
			xF	= abs(aux_xF);
			x01	= abs(x01);			
		}//End sentence if 
		
		
		//
		x0_prime = 1.0/(2.*(space_grid.k1[kmom] + av[0].f[kclock1[1]][1] ) );		
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<= 0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((x0_prime-x01)/ds)) ;
		
		
		
		//Loop on s variable 
		integral_on_s		= 0.0;
		double int_1		= 0.0;
		double int_2		= 0.0;	
		
		
		
		for (int kspace=0; kspace<1*Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x01 + kspace*ds;	
			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			taux1				= (g.t[kclock1[1]] + abs((s-x01)/(space_grid.k1[kmom])) - g.t[0]);
			aux_index			= floor(taux1/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			
			//Starting point 
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[1]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[1]][1];	
			//End sentece on time (time integration) 
			
			
			

			int_1+=  -  pow( soft_core + pow(integral_on_time + s, 2.), -1./2.)*ds/space_grid.k1[kmom] + 
						pow( soft_core + pow(s, 2. ), -1./2. )*ds/(space_grid.k1[kmom] + av[0].f[kclock1[1]][1] );
						
		}//End loop s variable
		
		
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>=-0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((xF-x0_prime)/ds)) ;
		
		
		
		
		
		//Loop on s variable 
		integral_on_s		= 0.0;
		for (int kspace=0; kspace<Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x0_prime + kspace*ds;			
			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			
			taux1				= (g.t[kclock1[1]] + abs((s-x0_prime)/space_grid.k1[kmom]) - g.t[0]);
			aux_index			= floor(taux1/dt);
			
			
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			
			//Starting sentence on time integration 
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[1]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[1]][1];	
			//End sentece on time (time integration) 
			
			
			int_2+= pow( soft_core + pow(integral_on_time + s, 2.), -1./2. )*ds/space_grid.k1[kmom];
						
		}//End loop s variable	
		
		
		integral_on_s	= int_1*1 + int_2;
		phase1[kmom]	= integral_on_s;
		
		
	}//End momentum loop
	
}//End Manfred's Formula integration on time*/




//Laser-Coulomb Coupling IR
void laser::laser_coulomb_coupling_x0II(grid &space_grid, double ds, double x01, double xF)
{
	//Here I integrate Manfred's formula for
	//Laser-Coulomb Coupling between IR laser and coulomb potential well
	//The formula is integrated on time 
	
	
	phase2 = (double*)malloc(space_grid.n1*sizeof(double));
	
	
	double s;
	double aux_xF		= xF;
	
	double integral_on_s	;
	double integral_on_time	;
	
	double taux1;
	double x0_prime		= 0.;
	
	double soft_core	= 0.;
	long double aux_index	;
	
	int Ns					;	
	int final_ktime			;	
	
	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {
		
		x0_prime = 1.0/(2.*( space_grid.k1[kmom] + av[0].f[kclock1[2]][1] ) );
		
		
		
		if (space_grid.k1[kmom]<0){
			ds	= -abs(ds);		
			x01	= -abs(x01);	
			xF	= -abs(aux_xF);	
		}
		else{
			ds	= +abs(ds);		
			x01	= abs(x01);		
			xF	= +abs(aux_xF);	
		}
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;
		else
			Ns   = floor(abs((x0_prime-x01)/ds)) ;
		
		
		
		//Loop on s variable 
		integral_on_s		= 0.0;
		double int_1		= 0.0;
		double int_2		= 0.0;	
		
		for (int kspace=0; kspace<1*Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x01 + kspace*ds;	
			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			taux1				= g.t[kclock1[2]] + abs((s-x01)/space_grid.k1[kmom]) - g.t[0];
			aux_index			= floor(taux1/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			
			//Starting point 
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[2]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[2]][1];	
			//End sentece on time (time integration) 
			
			
			
			int_1+=  - pow( soft_core + pow(integral_on_time+s,2.),-1./2.)*ds/space_grid.k1[kmom] + 
					 + pow( soft_core + pow(s,2.),-1./2.)*ds/( space_grid.k1[kmom] + av[0].f[kclock1[2]][1] );
			
		}//End loop s variable
		
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((xF-x0_prime)/ds)) ;
		
		
		
		//Loop on s variable  
		for (int kspace=0; kspace<Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x0_prime + kspace*ds;			
			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			
			taux1				= (g.t[kclock1[2]] + abs((s-x0_prime)/(space_grid.k1[kmom])) - g.t[0]);
			aux_index			= floor(taux1/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[2]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[2]][1];	
			//End sentece on time (time integration)
			
			
			
			int_2+= pow( soft_core + pow(integral_on_time+s,2.),-1./2.)*ds/space_grid.k1[kmom];
			
			
		}//End loop s variable	
		
		
		integral_on_s	= int_1*1 + int_2;
		phase2[kmom]	= integral_on_s;
		
		
	}//End momentum loop
	
}//End Manfred's Formula integration on time*/







//Laser-Coulomb Coupling IR 
void laser::laser_coulomb_coupling_fxI(  grid &space_grid, double ds, double x01, double xF )
{
	//Here I integrate Manfred's formula for
	//Laser-Coulomb Coupling between IR laser and coulomb potential well
	//The formula is integrated on time 
	
	
	phase1=(double*)malloc(space_grid.n1*sizeof(double));
	
	
	double s;
	double aux_xF=xF;
	double integral_on_s;
	double integral_on_time;
	double taux1;
	long double aux_index;
	double x_p_a = 0.;
	

	int Ns;	
	int Ns0;	
	int final_ktime;	
	x01		= abs(x01);

	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {		
		
		x_p_a	= 1.0/2./(space_grid.k1[kmom]+av_int[0].f[kclock1[1]][1]);		
		if (space_grid.k1[kmom]<0){
			ds	= -abs(ds);
			xF	= -abs(aux_xF);	
			x01 = -x01;
		}
		else{
			ds	= +abs(ds);
			xF	= +abs(aux_xF);
			x01 = +abs(x01);			
		}
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((x_p_a-x01)/ds)) ;
		
		
		
		//Loop on s variable 
		integral_on_s			= 0.;
		for (int kspace=0; kspace<Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x01 + kspace*ds;			
						
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			taux1				= (g.t[kclock1[1]] + abs((s-x01)/(space_grid.k1[kmom])) - g.t[0]);
			aux_index			= floor(taux1/dt);
						
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
						
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[1]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[1]][1];	
			//End sentece on time (time integration)
			

			
			
			//double soft_core = 0.5;
			integral_on_s+= -integral_on_time*ds/s/s; //
			
		}//End loop s variable
		phase1[kmom] = -integral_on_s/space_grid.k1[kmom];
		
	}//End momentum loop
	
}//End Manfred's Formula integration on time*/





//Laser-Coulomb Coupling IR 
void laser::laser_coulomb_coupling_fxII( grid &space_grid, double ds, double x02, double xF )
{
	//Here I integrate Manfred's formula for
	//Laser-Coulomb Coupling between IR laser and coulomb potential well
	//The formula is integrated on time 
	
	
	phase2=(double*)malloc(space_grid.n1*sizeof(double));
	
	
	double s;
	double aux_xF=xF;
	double integral_on_s;
	double integral_on_time;
	double taux1;
	long double aux_index;
	
	int Ns;	
	int final_ktime;	
	x02		=abs(x02);
	//Set_Av_Integral();	
	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {
		
			
		if (space_grid.k1[kmom]<0){
			ds = -abs(ds);
			xF = -abs(aux_xF);
			x02 = -x02;
		}
		else{
			ds = +abs(ds);
			xF = +abs(aux_xF);
			x02= abs(x02);
		}
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((xF-x02)/ds)) ;
		
		//Loop on s variable 
		integral_on_s			= 0.;
		for (int kspace=0; kspace<Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x02 + kspace*ds;			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			taux1				= (g.t[kclock1[2]] + abs((s-x02)/(space_grid.k1[kmom])) - g.t[0]);
			aux_index			= floor(taux1/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[2]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[2]][1];	
			//End sentece on time (time integration)
			
			
			//double soft_core = 0.5;
			integral_on_s+=  -integral_on_time*ds/s/s;
			
		}//End loop s variable
		phase2[kmom] = integral_on_s/abs(space_grid.k1[kmom]);
		
	}//End momentum loop
	
}//End Manfred's Formula integration on time*/




//Laser-Coulomb Coupling IR
void laser::laser_coulomb_coupling_x0(double *x0,grid &space_grid, double ds, double xF)
{
	//Here I integrate Manfred's formula for
	//Laser-Coulomb Coupling between IR laser and coulomb potential well
	//The formula is integrated on time 
	
	
	phase1=(double*)malloc(space_grid.n1*sizeof(double));
	phase2=(double*)malloc(space_grid.n1*sizeof(double));
	
	
	double x01;
	double s;
	double aux_xF=xF;
	double integral_on_s;
	double integral_on_s2;
	double integral_on_time;
	long double aux_index;
	
	int Ns;	
	int final_ktime;	
	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {
		
		//x01		= 0.1;//1./(2.*(space_grid.k1[kmom] + av[0].f[kclock1[1]][1] ) );		
		
		
		if (space_grid.k1[kmom]<0){
			ds = -abs(ds);
			xF = -abs(aux_xF);
			x01		= -0.01;
		}
		else{
			ds = +abs(ds);
			xF = +abs(aux_xF);
			x01		= 0.01;
		}
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;		
		else
			Ns   = floor(abs((xF-x0[kmom])/ds)) ;
		
		
		
		
		//Loop on s variable 
		integral_on_s = 0.;
		for (int kspace=0; kspace<Ns; kspace++) {			
			
			//Axis of the variable s
			s					= x0[kmom] + kspace*ds;			
			
			//Sentence on time (time integration)
			integral_on_time	= 0.;	
			final_ktime			= 0;
			
			double taux = (g.t[kclock1[1]] + abs((s-x0[kmom])/(space_grid.k1[kmom])) - g.t[0]);
			
			aux_index=floor(taux/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[1]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[1]][1];	
			//End sentece on time (time integration)
			
			
			//double soft_core = 0.5;
			integral_on_s+=  -integral_on_time*ds/s/s;
			//0*integral_on_time*integral_on_time*pow(s,-3.)*ds+ 
			//-0*integral_on_time*integral_on_time*integral_on_time*pow(s,-4.)*ds;
			
			
			
			
			
			//Second part or second phase using the pulse
			//x01 = 0.1;//1./(2.*(space_grid.k1[kmom] + av[0].f[kclock1[2]][1] ) );		
			
			//Loop on s variable kclock1[2]
			integral_on_s2		= 0.;			
			s					= x01 + kspace*ds;
			
			
			//Sentece on time
			integral_on_time	=0.;
			final_ktime			=0;
			aux_index			=floor((g.t[kclock1[2]] + 
										abs((s-x01)/(space_grid.k1[kmom])) - g.t[0])/dt);
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time= av_int[0].f[Nmaxt-1][1]-av_int[0].f[kclock1[2]][1];					
			else 
				integral_on_time= av_int[0].f[final_ktime-1][1]-av_int[0].f[kclock1[2]][1];
			//End sentece on time
			
			
			integral_on_s2+=  -integral_on_time*ds/s/s ;
			//-0*integral_on_time*integral_on_time*ds/s/s/s;
			
		}//End loop s variable
		phase1[kmom] = integral_on_s/abs(space_grid.k1[kmom]);
		phase2[kmom] = integral_on_s2/abs(space_grid.k1[kmom]);		
		
		
		
		if(kmom%(space_grid.n1/10)==0){
			cout << "\nIntegral1= " << integral_on_s;	
			cout << "     Integral2= " << integral_on_s2;			
		}
		
		
	}//End momentum loop
	
}//End Manfred's Formula integration on time*/





//Other Form for Induced-laser-Coupling
void laser::laser_coulomb_coupling_full(double xF, double ds, double soft_core, double Zatomic, grid &space_grid)
{
	
	phase1=(double*)malloc(space_grid.n1*sizeof(double));
	phase2=(double*)malloc(space_grid.n1*sizeof(double));
	
	
	double x01;
	double x02;
	double s;
	double xIR;
	double integral_on_s;
	double integral_on_s2;
	
	
	double integral_on_time;	
	double aux_xF		= xF;
	double tau_prime;
	long double aux_index;

	
	int Ns;
	int final_ktime; 	

	
	//Momentum loop
	for (int kmom=0; kmom<space_grid.n1; kmom++) {
		
		x01		= 1./(2.*(space_grid.k1[kmom] + av[0].f[kclock1[1]][1] ) );	
		x02		= 1./(2.*(space_grid.k1[kmom] + av[0].f[kclock1[2]][1] ) );			
		
		
		if (space_grid.k1[kmom]<0){
			ds = -abs(ds);
			xF = -abs(aux_xF);
		}
		else{
			ds = +abs(ds);
			xF = +abs(aux_xF);
		}
		
		
		
		//Number point on s variable
		if (space_grid.k1[kmom]>= -0.001 && space_grid.k1[kmom]<=0.001)
			Ns   = 150000;
		else
			Ns   = floor(abs((xF-x01)/ds)) ;
		
		
		
		
		
		//Loop on s variable 
		integral_on_s			= 0.;
		integral_on_s2			= 0.;		
		
		for (int kspace=0; kspace<Ns; kspace++) 
		{	
			
			integral_on_time	= 0.;
			s					= x01 + kspace*ds;			
			
						
			tau_prime			= abs(g.t[kclock1[1]] + abs((s-x01)/space_grid.k1[kmom]) - g.t[0]);			
			aux_index			= floor(tau_prime/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;
			else
				final_ktime=int(aux_index);
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[1]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[1]][1];	
			
						
			xIR			  = s + integral_on_time; 
			integral_on_s+= ds/sqrt(soft_core + xIR*xIR); 
			
			
			
			
			
			s			= x02 + kspace*ds;			
						
			tau_prime	= abs(g.t[kclock1[2]] + abs((s-x02)/space_grid.k1[kmom]) - g.t[0]);			
			aux_index	= floor(tau_prime/dt);
			
			
			if (aux_index>=2147483640)
				final_ktime=2147483641;			
			else
				final_ktime=int(aux_index);
			
			
			if (final_ktime>=Nmaxt-1)
				integral_on_time = av_int[0].f[Nmaxt-1][1] - av_int[0].f[kclock1[2]][1];
			else
				integral_on_time = av_int[0].f[final_ktime-1][1] - av_int[0].f[kclock1[2]][1];	
			
			
			xIR = s+integral_on_time; 
			integral_on_s2+= ds/sqrt(soft_core + xIR*xIR); 			
			
			
		}//End loop s variable
		
		phase1[kmom]= Zatomic*integral_on_s/abs(space_grid.k1[kmom]);
		phase2[kmom]= Zatomic*integral_on_s2/abs(space_grid.k1[kmom]);		
			
	}//End loop momentum

}//End function
	

	
	
//Finding index values 
int laser::vector_potential_match(double tI, double tII)
{
	//This function finds the index value correspondent to 
	//A(tI)=A(t?) NOT at tII A(tII)
	//tI the first time
	//tII is the upper limit to define the region to find A(t?)=A(tI)
	
	int second_ion_kont;
	int flag;
	int kstart = floor((tI-g.t[0])/dt);
	int kfinal = floor((tII-g.t[0])/dt);
	
	double Avn1  = av[0].f[kstart][1];
	double slope = (av[0].f[kstart+1][1]-av[0].f[kstart][1])/dt;
	
	
	//Choosing value for flag condition
	if (slope>0.)
		flag=0;
	else
		flag=1;
	
	
	
	//Selection of the parameters 
	for (int ktime=kstart; ktime<kfinal; ktime++)
	{	
	    slope= (av[0].f[ktime+1][1]-av[0].f[ktime][1])/dt;
		
		if (slope<0.) {
			if ((av[0].f[ktime][1]<Avn1)&(flag==0)) {
				second_ion_kont = ktime-1;
				flag            = 1;
			}			
		}
		else {
			if ((av[0].f[ktime][1]>Avn1)&(flag==1)) {
				second_ion_kont = ktime-1;
				flag            = 0;
			}			
		}	
	}//End loop for
	return second_ion_kont;
}//End finding index to the value


void laser::pulses_display()
{
	cout << "\n//*************************//\n";
	cout << "   LASER PARAMETERS ";	
	cout << "\n//*************************//\n";


	for (int kpulses=0; kpulses<Npulses; kpulses++) {
		
		cout << "\n\n//******  Laser Parameters for pulse Nº" << kpulses  <<"  ******//";
		cout << "\nIntensity"	<< kpulses << "= "	  <<	I0[kpulses] ;
		cout << "\nCFrequency"	<< kpulses << "= "    <<	w0[kpulses] ;
		cout << "\nCycles"		<< kpulses << "= "	  <<	cycles0[kpulses];
			cout << "   BaseWidth"   << kpulses << "= "    <<    twidth[kpulses];
		cout << "\nChirp"		<< kpulses << "= "	  <<	chirp[kpulses];
		cout << "\nEllipticity" << kpulses << "= "	  <<	e[kpulses];
		cout << "\nEnvelope"	<< kpulses << "= "	  <<	envelope[kpulses];

		if (kpulses>0) 
			cout << "\n\nDelay"	<< kpulses-1 << "= "	  <<	delay0[kpulses-1];
		cout << "\n\nTimeAtStart" << kpulses << "= "	  <<	clock0[kpulses];
		cout << "\nTimeAtMaxima"<< kpulses << "= "	  <<	clock1[kpulses];
		cout << "\nTimeAtEnd"   << kpulses << "= "	  <<	clock2[kpulses];				
	}
	cout << "\n\n//******* End laser parameters *******//\n\n";
}



//Saving function for laser pulses: electric field and vector potential
void laser::saving_pulses(FILE *timeaxis,FILE *file, FILE *integrals)
{

	for (int ktime=0; ktime<Nmaxt; ktime++)
		fprintf(timeaxis,"%e \n",g.t[ktime]);
	
	
	for (int ktime=0; ktime<Nmaxt; ktime++) 
		fprintf(file,"%e \n",efield.f[ktime][1]);
	
	
	for (int ktime=0; ktime<Nmaxt; ktime++) 
		fprintf(file,"%e \n",avector.f[ktime][1]);	


	//Saving electric field pulses
	for (int kpulses=0; kpulses<Npulses; kpulses++) 
		for (int ktime=0; ktime<Nmaxt; ktime++)
			fprintf(file,"%e \n",ef[kpulses].f[ktime][1]);
		

	//Saving electric field envelopes	
	for (int kpulses=0; kpulses<Npulses; kpulses++) 
		for (int ktime=0; ktime<Nmaxt; ktime++)
					fprintf(file,"%e \n",env[kpulses].f[ktime][1]);
	
	
	
	//Saving vector potentials
	for (int kpulses=0; kpulses<Npulses; kpulses++) 
		for (int ktime=0; ktime<Nmaxt; ktime++)
			fprintf(file,"%e \n",av[kpulses].f[ktime][1]);
			
	
	
	//Saving Integral of the vector potential
	//Set_Av_Integral();
	for (int kpulses=0; kpulses<Npulses; kpulses++) 
		for (int ktime=0; ktime<Nmaxt; ktime++) 
			fprintf(integrals,"%e %e\n",av_int[kpulses].f[ktime][1],
					avsq_int[kpulses].f[ktime][1]);
		
	
}//End saving function


void laser::save_azeros(FILE *file)
{
	for (int kzeros=0; kzeros<2*cycles0[0]+1; kzeros++)
		fprintf(file,"%d \n",Azeros[kzeros]);
}


void laser::save_ezeros(FILE *file)
{
	for (int kzeros=0; kzeros<2*cycles0[0]+2; kzeros++)
		fprintf(file,"%d \n",Ezeros[kzeros]);
}


/*========================== END  SECUNDARY FUNCTIONS =============================*/

#endif
//END
