#ifndef DIAGONALIZE1D_H
#define DIAGONALIZE1D_H
#include <iostream>
//using namespace std;
class diagonalize1D: public fullhamiltonian
{	
public:	
	
	
	char jobz;
	char range;			
	
	
	
	MKL_INT ldz;
	MKL_INT matrix_order;
	
	
	
	
	
	MKL_INT il;		 //Controling initial index number of eigen value and eigen-vector to calculate
	MKL_INT iu;		//Controling final index number of eigen value and eigen-vector to calculate
	MKL_INT nzc;	
	
	
	
	
	MKL_INT *isuppz;
	MKL_INT *m;		
	MKL_INT *tryrac;

	
	
	
	double vl;
	double vu;
	
		
	
	void initialize(){		

		ldz				= Nx;
		matrix_order	= LAPACK_COL_MAJOR;		//101;//
		
		
		jobz			= 'V';
		range			= 'I';	
		
		
		vl				= .1;
		vu				= .8;
		
		il				= 1;				//Controling initial index number of eigen value and eigen-vector to calculate
		iu				= MKL_INT(Neig);    //Controling final index number of eigen value and eigen-vector to calculate		
		nzc				= MKL_INT(Neig);
		
		isuppz			= (MKL_INT*)mkl_malloc(2*Neig*sizeof(MKL_INT),16);		//new MKL_INT[2*neig];
		m				= (MKL_INT*)mkl_malloc((iu-il+1)*sizeof(MKL_INT),16);	//new MKL_INT[iu-il+1];
		tryrac			= (MKL_INT*)mkl_malloc(Neig*sizeof(MKL_INT),16);		//new lapack_logical[Neig];								//(MKL_INT*)mkl_malloc(Neig*sizeof(MKL_INT),16)  ;		//new MKL_INT[neig];		
		
		
		memset(isuppz, 0., 2*Neig );
		memset(m,  0., (iu-il+1) );
		memset(tryrac,  0., Neig );		
		
		//cout << "\nNx="<<Nx;
	}
	

	
	//Building eigen vectors 
	void diagonalization(wavefunction *w){			
		
		vectors_diagonals();
		MKL_INT int0	= LAPACKE_dstemr(matrix_order, jobz, range, Nx, diag, up_diag, vl, vu, il, iu, m, e_values, e_vectors, ldz, nzc, isuppz, tryrac);
				
		
		//Writing data on wavefunction
		for (MKL_INT j=0; j<Neig; j++) {
			
			for (MKL_INT i=0; i<Nx; i++) {
				
				w[j].w[i][0] = e_vectors[index12_2d(j,i)];
				w[j].w[i][1] = 0.;
				
			}			
			
		}//*/

	}
	

	void diagonalization_length_g(wavefunction *w, double _efield1){
		
		
		//Bulding symmetrics vectors 
		vectors_diagonals_dynamic(_efield1);

		//Solving eigen values and eigen vectors problem 
		MKL_INT int0	= LAPACKE_dstemr(matrix_order, jobz, range, Nx, diag, up_diag, vl, vu, il, iu, m, e_values, e_vectors, ldz, nzc, isuppz, tryrac);
		
		
		
		//Loop for Writing data on wavefunction
		for (MKL_INT j=0; j<Neig; j++) {
			
			for (MKL_INT i=0; i<Nx; i++) {
				
				w[j].w[i][0] = e_vectors[index12_2d(j,i)];
				w[j].w[i][1] = 0.;
				
			}
			
		}//End loop*/
		
	}//End */
	
	
	
	// Destructor	
	~diagonalize1D()
	{ 
		mkl_free(isuppz);
		mkl_free(m);
		mkl_free(tryrac);
	}	//*/		

	

};
#endif
//END 
