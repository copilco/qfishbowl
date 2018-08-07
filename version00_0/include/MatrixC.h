  /*
  *
  *    Created by Rafael Moran on 24/05/2009.
  *
  */


#ifndef MATRIX_C
#define MATRIX_C

//Start Class that implement a complex number
class NC{ 
	public:
		double r;
		double i;
};//End Class complex number


//Start Class that implement a complex matrix
class MatrixC{
	    int nf;	
	    int nc;	
	    NC * vector;

	public:   
	    MatrixC(int nf, int nc);
	    ~MatrixC();             

	    NC Get(int i, int j);
	    void Set(int i, int j, double vr, double vi);

};//End Class complex matrix


//Start of variables
MatrixC::MatrixC(int _nf, int _nc) 
{			   

    vector = new NC[_nf*_nc]; 
    nf = _nf;
    nc = _nc;

    for (int i=0; i<nf*nc; i++)
	{
        vector[i].r = 0.00; 
        vector[i].i = 0.00; 
    }
}//End of start variables


//Destructor Matrix
MatrixC::~MatrixC(){
    delete vector;
}


//Get values of vector
NC MatrixC::Get(int i, int j){
    return vector[i*nc+j];
}


//Start function Set
void MatrixC::Set(int i, int j, double vr, double vi)
{
    	vector[i*nc+j].r = vr;
	vector[i*nc+j].i = vi;
}//End Set

#endif
//End
