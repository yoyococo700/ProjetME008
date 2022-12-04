#define n 100


/*
    d   a   0   0
C = a   d   a   0
    0   a   d   a
    0   0   a   d


    l1  0   0   0
L = u2  l2  0   0
    0   u3  l3  0
    0   0   u4  l4

*/




//Matrices de Toeplitz tridiagonales
void factorisation_Cholesky(float a,float d, float* l,float* u){
    
}





void resol_trig_inf( float **A, float x[], float b[]) {
	int i,k;
	float l;
	
	for (i=0;i<n;i++) {
		l=0;
		for (k=0;k<=i-1;k++) l += A[i][k]*x[k];
		x[i] = (b[i] -l)/A[i][i];
	}
}


void resol_trig_sup( float **A , float x[], float b[]) {
	int i,k;
	float l;
	
	for (i=n;i>=0;i--) {
		l=0;
		for (k=i+1;k<n;k++) l += A[i][k]*x[k];
		x[i] = (b[i] -l)/A[i][i];
	}
}	


