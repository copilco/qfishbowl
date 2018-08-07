/*
 *  laser.cpp laser pulses train differents characteristics by pulse
 *
 *  Created by Alexis Chacón laser0_0_1
 *
 */

#ifndef LASER0_0_1_H
#define LASER0_0_1_H
#include <stdlib.h>
#include <string>
#include "timegrid.h"
#include "timeobject.h"
#include "grid.h"
#include "constant.h"
#include <vector>
using namespace std;

class laser0_0_1 
{

public:
	
	char *check;
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
	
	
	
	/*==========================*/
	/*      MAIN FUNCTIONS
	/*==========================*/
	laser0_0_1(int _Npulses);      							                                // Creator Object
	//~laser();																			//Destructor
	void laser_pulses(double _dt, double _t01, double _blaser, double _alaser);	        // LASER PULSES TRAIN
	
	/*=========================*/
	/*   SECUNDARY FUNCTIONS
	/*=========================*/
	inline void initialize_amplitude_period();			// Initialize max amplitude and period
	inline void set_start_time_end_time();		        // "Start" and "end" time per pulse
	void laser_time_grid();									// Object time axis and field_set
	void evaluating_laser_pulse();						// Evaluation pulses train
	void adding_laser_pulses();									// Pulses sum
	void set_vector_potential();						// Vector potential per pulse
	void vector_potential();							// Total vector potential
	void put_on_envelope();								// Generator of Envelope of the pulses
	void set_vector_potential_integral();
	inline void set_main_indexs_pulses();	
	void checking_starting_axis_time();
	
	double Lqmajor(vector<double>& v);					//Finding the maximum value of a vector
	double Lqminus(vector<double>& v);					//Finding the minimum value ot a vector
	
	
	
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
laser0_0_1::laser0_0_1(int _Npulses)
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
void laser0_0_1::laser_pulses(double _dt, double _t01, double _blaser, double _alaser)
{
	//return;
   	/*======= Pulse's Parameter ======*/
   	t01		= _t01;										//	Start time first pulse
   	dt		= abs(_dt);			                        //	Time step 	
   	
	blaser	= abs(_blaser);								//	Time before the laser or the pulse train
   	alaser	= abs(_alaser);								//	Time after the laser or the pulse train 	
	
	
   	//================================//
	initialize_amplitude_period();						//	Initialize max amplitude and period        
   	set_start_time_end_time();							//	"Start" and "end" time by each pulse
	
	
   	major0	= Lqmajor(clock2);	     	                //	Maximum time of all pulse
   	minus0	= Lqminus(clock0);							//	Minimun time of all pulse 
	
		
	//CHECKING CONTROL OF AXIS 
	checking_starting_axis_time();	
	set_main_indexs_pulses();		
	
	
	
   	laser_time_grid();										//	Objet time axis      
   	evaluating_laser_pulse();							//	Evaluation pulses
	
   	
	
	adding_laser_pulses();										//	Pulses sum (build pulses train)
	set_vector_potential();								//	Vector potential by each pulse
   	vector_potential();									//	Vector potential to pulses train 
	

	

}//End laser_pulses

/*============ END MAIN FUNCTIONS ===============*/



/*=================================================================*/
		/*========= SECUNDARY FUNCTIONS ===========*/

//== Function initialize max amplitude and period ==// 
inline void laser0_0_1::initialize_amplitude_period()
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
inline void laser0_0_1::set_start_time_end_time()
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


inline void laser0_0_1::set_main_indexs_pulses()
{
	double min_aux0 = 0;

	for (int kpulse=0;kpulse<Npulses;kpulse++)
	{		
		kclock0[kpulse] = floor( abs((clock0[kpulse]) - (minus0-blaser))/dt)	+ 1;
		kclock1[kpulse] = floor( abs((clock1[kpulse]) - (minus0-blaser))/dt )	+ 1;
		kclock2[kpulse] = floor( abs((clock2[kpulse]) - (minus0-blaser))/dt )	+ 1;			
	}//End loop

}

void laser0_0_1::checking_starting_axis_time(){
	
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
void laser0_0_1::laser_time_grid()   
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
	cout << "        About laser0_0_1.h ";	
	cout << "\nNtime= "<<Nmaxt;
	cout << "      Memory= "<<((2+3*Npulses)*24+8)*Nmaxt*1e-6 << "  Mb"<<endl;
	cout << "/********************************/\n\n\n";		
}


//== Funtion that evaluate of the pulses train ==//
void laser0_0_1::evaluating_laser_pulse()
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



void laser0_0_1::put_on_envelope(){
	
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
void laser0_0_1::adding_laser_pulses()
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
void laser0_0_1::set_vector_potential()
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
void laser0_0_1::set_vector_potential_integral()
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
void laser0_0_1::vector_potential()
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
double laser0_0_1::Lqmajor(vector<double>& v)
{
	int n=v.size();
	double may = v[0];
	
	for (int k=1;k<n;k++)   
		if (v[k]>may) may=v[k];
	
	return may;
}//End major



//Start function minus
double laser0_0_1::Lqminus(vector<double>& v)
{
	int n=v.size();
	double min=v[0];
	
	for (int k=1;k<n;k++)
		if (v[k]<min) min=v[k];	
	
	return min;
}//End minus






//Member function to calculate zeros on the electric field
void laser0_0_1::zeros_electric_f()
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
void laser0_0_1::zeros_vector_pot()
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



	
	
//Finding index values 
int laser0_0_1::vector_potential_match(double tI, double tII)
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




void laser0_0_1::pulses_display()
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
void laser0_0_1::saving_pulses(FILE *timeaxis,FILE *file, FILE *integrals)
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
	//set_vector_potential_integral();
	for (int kpulses=0; kpulses<Npulses; kpulses++) 
		for (int ktime=0; ktime<Nmaxt; ktime++) 
			fprintf(integrals,"%e %e\n",av_int[kpulses].f[ktime][1],
					avsq_int[kpulses].f[ktime][1]);
		
	
}//End saving function


void laser0_0_1::save_azeros(FILE *file)
{
	for (int kzeros=0; kzeros<2*cycles0[0]+1; kzeros++)
		fprintf(file,"%d \n",Azeros[kzeros]);
}


void laser0_0_1::save_ezeros(FILE *file)
{
	for (int kzeros=0; kzeros<2*cycles0[0]+2; kzeros++)
		fprintf(file,"%d \n",Ezeros[kzeros]);
}


/*========================== END  SECUNDARY FUNCTIONS =============================*/

#endif
//END
