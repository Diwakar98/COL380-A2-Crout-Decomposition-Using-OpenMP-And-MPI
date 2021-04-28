#include <omp.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>

int T;

void write_output(char fname[], double** arr, int n ){
	FILE *f = fopen(fname, "w");
	for( int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			fprintf(f, "%0.12f ", arr[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void crout(double **A, double **L, double **U, int n) {
	int i, j, k;
	double sum = 0;
	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}
}	

void crout1(double **A, double **L, double **U, int n) {
    #pragma omp parallel for num_threads(T)
	for (int i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	for (int j = 0; j < n; j++) {
        #pragma omp parallel for 
        for (int i = j; i < n; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
        #pragma omp parallel for 
		for (int i = j; i < n; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}
}	

void crout2(double **A, double **L, double **U, int n) {
	#pragma omp parallel sections num_threads(T)
	{
		#pragma omp section
		for (int i = 0; i < n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = n/16; i < 2*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 2*n/16; i < 3*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 3*n/16; i < 4*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 4*n/16; i < 5*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 5*n/16; i < 6*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 6*n/16; i < 7*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 7*n/16; i < 8*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 8*n/16; i < 9*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 9*n/16; i < 10*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 10*n/16; i < 11*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 11*n/16; i < 12*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 12*n/16; i < 13*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 13*n/16; i < 14*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 14*n/16; i < 15*n/16; i++) {
			U[i][i] = 1;
		}
		#pragma omp section
		for (int i = 15*n/16; i < n; i++) {
			U[i][i] = 1;
		}
	}
	for (int j = 0; j < n; j++) {
		double sum1 = 0;
		for (int k = 0; k < j; k++) {
			sum1 = sum1 + L[j][k] * U[k][j];	
		}
		L[j][j] = A[j][j] - sum1;
		U[j][j] = (A[j][j] - sum1) / L[j][j];
		
		#pragma omp parallel sections 
		{
		#pragma omp section
        for (int i = j+1 + 0*(n-j-1)/8; i < j+1 + 1*(n-j-1)/8; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp section
        for (int i = j+1 + 1*(n-j-1)/8; i < j+1 + 2*(n-j-1)/8; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp section
        for (int i = j+1 + 2*(n-j-1)/8; i < j+1 + 3*(n-j-1)/8; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp section
        for (int i = j+1 + 3*(n-j-1)/8; i < j+1 + 4*(n-j-1)/8; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp section
        for (int i = j+1 + 4*(n-j-1)/8; i < j+1 + 5*(n-j-1)/8; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp section
        for (int i = j+1 + 5*(n-j-1)/8; i < j+1 + 6*(n-j-1)/8; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp section
        for (int i = j+1 + 6*(n-j-1)/8; i < j+1 + 7*(n-j-1)/8; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp section
        for (int i = j+1 + 7*(n-j-1)/8; i < j+1 + 8*(n-j-1)/8; i++) {
			double sum = 0;
			for (int k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp section
        for (int i = j+1 + 0*(n-j-1)/8; i < j+1 + 1*(n-j-1)/8; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
		#pragma omp section
        for (int i = j+1 + 1*(n-j-1)/8; i < j+1 + 2*(n-j-1)/8; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
		#pragma omp section
        for (int i = j+1 + 2*(n-j-1)/8; i < j+1 + 3*(n-j-1)/8; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
		#pragma omp section
        for (int i = j+1 + 3*(n-j-1)/8; i < j+1 + 4*(n-j-1)/8; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
		#pragma omp section
        for (int i = j+1 + 4*(n-j-1)/8; i < j+1 + 5*(n-j-1)/8; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
		#pragma omp section
        for (int i = j+1 + 5*(n-j-1)/8; i < j+1 + 6*(n-j-1)/8; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
		#pragma omp section
        for (int i = j+1 + 6*(n-j-1)/8; i < j+1 + 7*(n-j-1)/8; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
		#pragma omp section
        for (int i = j+1 + 7*(n-j-1)/8; i < j+1 + 8*(n-j-1)/8; i++) {
			double sum = 0;
			for(int k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				exit(0);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
		
		}
	}
}	

void crout3(double **A, double **L, double **U, int n) {
	#pragma omp parallel for num_threads(T)
	for (int i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	for (int j = 0; j < n; j++) {
        double sum1 = 0;
		for (int k = 0; k < j; k++) {
			sum1 = sum1 + L[j][k] * U[k][j];	
		}
		L[j][j] = A[j][j] - sum1;
		U[j][j] = (A[j][j] - sum1) / L[j][j];
		#pragma omp parallel sections
		{
			#pragma omp section
			{
			#pragma omp parallel for 
			for (int i = j+1; i < n; i++) {
				double sum = 0;
				for (int k = 0; k < j; k++) {
					sum = sum + L[i][k] * U[k][j];	
				}
				L[i][j] = A[i][j] - sum;
			}
			}
			#pragma omp section
			{
			#pragma omp parallel for 
			for (int i = j+1; i < n; i++) {
				double sum = 0;
				for(int k = 0; k < j; k++) {
					sum = sum + L[j][k] * U[k][i];
				}
				if (L[j][j] == 0) {				
					exit(0);
				}
				U[j][i] = (A[j][i] - sum) / L[j][j];
			}
			}
		}
	}
}	

int main(int argc, char* argv[]){ 
    if (argc<4){
        printf("ERROR\n");
        return 0;
    }
    int n = atoi(argv[1]);
    char *filename = argv[2];
    T = atoi(argv[3]);
    int strat = atoi(argv[4]);
    // omp_set_num_threads(T);
	omp_set_nested(1);
    FILE *fptr = fopen(filename,"r");
    double **A = (double **)malloc(n*sizeof(double *));
    for(int i=0;i<n;i++){
        A[i] = (double *)malloc(n*sizeof(double));
    }
    double **L = (double **)malloc(n*sizeof(double *));
    for(int i=0;i<n;i++){
        L[i] = (double *)malloc(n*sizeof(double));
    }
    double **U = (double **)malloc(n*sizeof(double *));
    for(int i=0;i<n;i++){
        U[i] = (double *)malloc(n*sizeof(double));
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            fscanf(fptr,"%lf",&A[i][j]);
        }
    }
    fclose(fptr);
    if(strat==0){
        crout(A,L,U,n);
	}
    else if(strat==1){
        crout1(A,L,U,n);
	}
	else if(strat==2){
        crout2(A,L,U,n);
	}
	else if(strat==3){
        crout3(A,L,U,n);
	}
	char linit[] = "output_L_";
	char uinit[] = "output_U_";
	char score[] = "_";
	char ext[] = ".txt";
	char* lfile = (char *)malloc(300);
	lfile[0] = '\0';
	strcat(lfile,linit); strcat(lfile,argv[4]); strcat(lfile,score); strcat(lfile,argv[3]); strcat(lfile,ext);
	char* ufile = (char *)malloc(300);
	ufile[0] = '\0';
	strcat(ufile,uinit); strcat(ufile,argv[4]); strcat(ufile,score); strcat(ufile,argv[3]); strcat(ufile,ext);
	write_output(lfile,L,n);
	write_output(ufile,U,n);
	return 0;
}