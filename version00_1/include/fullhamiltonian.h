#ifndef FULLHAMILTONIAN_H
#define FULLHAMILTONIAN_H

#include <iostream>
#include "mkl.h"
#include "mkl_lapacke.h"
#include "grid.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "constant.h"
#include <assert.h>

#ifndef lapack_int
#define lapack_int MKL_INT
#endif

#ifndef lapack_logical
#define lapack_logical lapack_int
#endif


using namespace std;
class fullhamiltonian: public hamiltonian
{
	public:
	
	MKL_INT Nx;
	MKL_INT Neig;	

	double *diag;
	double *up_diag;
	
	double *e_vectors;	
	double *e_values;	
	
	
	//Iniatilizing arrays and allocating memory
	void initialize(){
		
		
		diag			=  (double*)mkl_malloc(Nx*sizeof(double),16);	
		up_diag			=  (double*)mkl_malloc((Nx-1)*sizeof(double),16);
			
		e_vectors		=  (double*)mkl_malloc(Neig*Nx*sizeof(double),16);	
		e_values		=  (double*)mkl_malloc(Neig*sizeof(double),16);		
		
		
		memset(e_vectors, 0., Neig*Nx );
		memset(e_values,  0., Neig );	
		
		/*for (int i=0; i<Neig*Nx; i++)
			e_vectors[i]=0.;
		
		for (int i=0; i<Neig; i++)
			e_values[i]=0.;

	//*/
	}
	
	
	//Loading variables
	void put_on_variables(int _Neig)
	{
		Nx		= MKL_INT( n1);		
		Neig	= MKL_INT(_Neig);	
		
		initialize();
	}
	
	// Building diagonals vectors: main and down_up diagonals	
	void vectors_diagonals()
	{
		
		for (MKL_INT i=0;i<Nx;i++)
		{

			diag[i]= 1.0/dx1/dx1 + v[i];
			
			if (i<Nx-1)
				up_diag[i]= -1.0/dx1/dx1/2.0;		
			
		}
	
	}//End vectors_diagonals

	
	// Building diagonals dynamic vectors: main and down_up diagonals	
	void vectors_diagonals_dynamic(double efield1)
	{
		
		for (MKL_INT i=0;i<Nx;i++)
		{
			
			diag[i]= 1.0/dx1/dx1 + v[i]+ efield1*x1[i];//+ efield1*h.x1[i];
			
			if (i<Nx-1)
				up_diag[i]= -1.0/dx1/dx1/2.0;		
			
		}
		
	}//End vectors_diagonals	
	
	
	// Destructor	
	~fullhamiltonian()
	{ 
		mkl_free(diag);
		mkl_free(up_diag);
		mkl_free(e_values);
		mkl_free(e_vectors);
				
	}	//*/
	


};
#endif
//END 
