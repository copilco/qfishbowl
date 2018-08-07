/*** Library Ordinary Differential Equations Solution Euler, Best Euler and Runge Kutta Fourth order method.
     This last method is implemented in zero, one and two dimention ----------------------------
     ....... Chacón, Alexis nov, 2009 .......***/
#ifndef DEQUATION_H
#define DEQUATION_H
#include <complex>
#define complex complex<double>
#include <vector>
#include <math.h>

using namespace std;

class dequation: public timegrid
{
public:
	int Nord;
	int n1;
	int n2;
	
	int a1;
	int a2;
	int ktime1;
	
	double h;   
	double hh;
	double a;
	double b;

	vector <complex> y;
	vector <complex> y1;
	vector <complex> y2;
	
	vector <complex> alpha;
	vector <complex> alpha1;
	vector <complex> alpha2;

	vector<complex> aux;
	vector <complex> kalpha0;
	vector <complex> kalpha1;
	vector <complex> kalpha2;
	vector <complex> kalpha3;
	
	dequation(int _Nord, int _n1, int _n2);
	void put_on_grid(timegrid &g);

	void Euler(complex (*function)(double, complex ) ,int ktime);
	void EulerB(complex (*function)(double, complex ) ,int ktime);
	void RungeK4( complex (*function)(double, complex ), int ktime );
	void RungeK4_1D( complex(*function)(int, int, double, complex ), int ktime );  
	void RungeK4_2D( complex(*function)(int, int, int, double, complex ), int ktime );
	void RungeK4C_2D( void (*function)( int, double, vector<complex>& ), int ktime );
	
	int index2(int p1, int p2);
	int index3(int p1, int p2, int p3);
};


dequation::dequation(int _Nord, int _n1, int _n2)
{
	Nord=_Nord;
	n1=_n1;   
	n2=_n2;
	
	alpha.resize(Nord, 0 );
	alpha1.resize(Nord*n1, 0 );
	alpha2.resize(Nord*n1*n2, 0 );

	aux.resize(n1*n2,0);
	kalpha0.resize(n1*n2,0);
	kalpha1.resize(n1*n2,0);
	kalpha2.resize(n1*n2,0);
	kalpha3.resize(n1*n2,0);
}


void dequation::put_on_grid(timegrid &g)
{
	set_grid(g.n, g.dt, g.t0);
	symmetry_n=g.symmetry_n;
	
	y.resize(n, 0 );
	y1.resize(n1, complex(0,0));
	y2.resize(n1*n2, complex(0,0));
	
	h  = g.dt;
	hh = 0.5*g.dt;
	a  = 2.0;
	b  = 1.0/6.0;
}


void dequation::Euler(complex (*function)(double, complex ) ,int ktime){
	alpha[0] = (*function)( t[ktime], y[ktime] );
 	y[ktime+1]=y[ktime]+h*alpha[0];
}


void dequation::EulerB(complex (*function)(double, complex ), int ktime ){
	alpha[0] = (*function)( t[ktime], y[ktime] );
	alpha[1] = (*function)( (t[ktime]+h), (y[ktime]+h*alpha[0]) );
 	y[ktime+1]=y[ktime]+hh*(alpha[0]+alpha[1]);
}


void dequation::RungeK4( complex (*function)(double, complex ), int ktime)
{    	                
	alpha[0] = h*(*function)( t[ktime], y[ktime] );
	alpha[1] = h*(*function)( (t[ktime]+hh), ( y[ktime] + alpha[0]/a ) );
	alpha[2] = h*(*function)( (t[ktime]+hh), ( y[ktime] + alpha[1]/a ) );
	alpha[3] = h*(*function)( (t[ktime]+ h), ( y[ktime] + alpha[2]) );  
	
	y[ktime+1] = y[ktime] + b*( alpha[0] +  a*( alpha[1]+ alpha[2] ) + alpha[3] );
}


void dequation::RungeK4_1D( complex(*function)(int, int, double, complex ), int ktime )
{
	ktime1=2*ktime;
	
	for (int i=0;i<n1;i++)
		alpha1[index2(0,i)] = h*(*function)( i, ktime1, t[ktime], y1[i] );
	
	ktime1+=1;
	
 	for (int i=0;i<n1;i++)
		alpha1[index2(1,i)] = h*(*function)( i, ktime1, (t[ktime]+hh), (y1[i] + alpha1[index2(0,i)]/a) );

	for (int i=0;i<n1;i++)
		alpha1[index2(2,i)] = h*(*function)( i, ktime1, (t[ktime]+hh), (y1[i] + alpha1[index2(1,i)]/a) );
	
	ktime1+=1;
	
	for (int i=0;i<n1;i++)
		alpha1[index2(3,i)] = h*(*function)( i, ktime1, (t[ktime]+ h), (y1[i] + alpha1[index2(2,i)]) );  
        
	for (int i=0;i<n1;i++)
		y1[i] +=b*( alpha1[index2(0,i)]  +  a*( alpha1[index2(1,i)] +
											   alpha1[index2(2,i)] ) + alpha1[index2(3,i)] ); 
}


void dequation::RungeK4_2D( complex(*function)(int, int, int, double, complex ), int ktime )
{
	ktime1=2*ktime;
	
	for (int j=0;j<n2;j++)
		for (int i=0;i<n1;i++)
			alpha2[index3(0,j,i)] = h*(*function)(i, j, ktime1, t[ktime], y2[index2(j,i)] );          
	
	ktime1+=1;
	
	for (int j=0;j<n2;j++)
		for (int i=0;i<n1;i++)
			alpha2[index3(1,j,i)] = h*(*function)(i, j, ktime1, (t[ktime]+hh), (y2[index2(j,i)] + alpha2[index3(0,j,i)]/a ) );          
	
	for (int j=0;j<n2;j++)
		for (int i=0;i<n1;i++)
	        alpha2[index3(2,j,i)] = h*(*function)(i, j, ktime1, (t[ktime]+hh), (y2[index2(j,i)] + alpha2[index3(1,j,i)]/a ) );          
	
	ktime1+=1;
	
	for (int j=0;j<n2;j++)
		for (int i=0;i<n1;i++)
			alpha2[index3(3,j,i)] = h*(*function)(i, j, ktime1, (t[ktime]+ h), (y2[index2(j,i)] + alpha2[index3(2,j,i)]) );            
	
	for (int j=0;j<n2;j++)
		for (int i=0;i<n1;i++)
			y2[index2(j,i)]  += b*( alpha2[index3(0,j,i)] + a*( alpha2[index3(1,j,i)] +
															   alpha2[index3(2,j,i)] ) + alpha2[index3(3,j,i)] );	    
}


void dequation::RungeK4C_2D( void (*function)( int, double, vector<complex>& ), int ktime )
{   
	ktime1=2*ktime;
	(*function)(ktime1, t[ktime], y2);

	for (int j=0;j<n1*n2;j++){
		kalpha0[j] = h*kalpha3[j];
		aux[j] = y2[j] + kalpha0[j]/a;
	}
	
	ktime1+=1;
	(*function)(ktime1, t[ktime]+hh, aux);
	
	for (int j=0;j<n1*n2;j++){
		kalpha1[j] = h*kalpha3[j];
		aux[j] = y2[j] + kalpha1[j]/a;	
	}
	
	(*function)(ktime1, t[ktime]+hh, aux);
	
	for (int j=0;j<n1*n2;j++){
		kalpha2[j] = h*kalpha3[j];
		aux[j] = y2[j] + kalpha2[j]; 	
	}
	
	ktime1+=1;
	(*function)(ktime1, t[ktime]+h, aux );
	
	for (int j=0;j<n1*n2;j++)
		y2[j] +=  b*( kalpha0[j] + a*( kalpha1[j] + kalpha2[j] ) + h*kalpha3[j] );	    
}


int dequation::index2(int p1, int p2)
{
	a1 = n1*p1+p2;
	return a1;  
}


int dequation::index3(int p1, int p2, int p3)
{
	a2 = n1*n2*p1 + index2(p2,p3);//n1*p2+p3;
	return a2;  
}
#endif
