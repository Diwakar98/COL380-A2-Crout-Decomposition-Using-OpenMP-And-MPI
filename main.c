#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>

int number_of_threads;

/*
void crout0(double **A, double **L, double **U, int n) {
    int i, j, k;
    double sum = 0;
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (j = 0; j < n; j++) {
        printf("j : %d\n",j);
        sum = 0;
        // printf("\ti : %d\n",j);
        printf("------------------------------\n");
        for (k = 0; k < j; k++) {
            sum = sum + L[j][k] * U[k][j];
            // printf("\t\tREAD: \tL[%d,%d]\n",j,k);
            // printf("\t\tREAD: \tU[%d,%d]\n",k,j);
        }
        L[j][j] = A[j][j] - sum;
        for (i = j+1; i < n; i++) {
            sum = 0;
            for (k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];
                printf("\t\tREAD: \tL[%d,%d]\n",i,k);
                printf("\t\tREAD: \tU[%d,%d]\n",k,j);
            }
            L[i][j] = A[i][j] - sum;
            printf("\t\tWRITE: \tL[%d,%d]\n",i,j);
        }
        printf("--------\n");
        for (i = j; i < n; i++) {
            sum = 0;
            // printf("\ti : %d\n",i);
            for(k = 0; k < j; k++) {
                sum = sum + L[j][k] * U[k][i];
                printf("\t\tREAD: \tL[%d,%d]\n",j,k);
                printf("\t\tREAD: \tU[%d,%d]\n",k,i);
            }
            if (L[j][j] == 0) {				
                exit(0);
            }
            U[j][i] = (A[j][i] - sum) / L[j][j];
            printf("\t\tREAD: \tL[%d,%d]\n",j,j);
            printf("\t\tWRITE: \tU[%d,%d]\n",j,i);
        }
    }
}	
*/

void crout0(double **A, double **L, double **U, int n) {
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
/*
void crout1(double **A, double **L, double **U, int n, int number_of_threads) {
    int i, j, k;
    double sum = 0;
    #pragma parallel for num_threads(number_of_threads)
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (j = 0; j < n; j++) {
        #pragma parallel for num_threads(number_of_threads)
        for (i = j; i < n; i++) {
            printf("X\n");
            sum = 0;
            for (k = 0; k < j; k++) {
                sum = sum + L[i][k] * U[k][j];	
            }
            L[i][j] = A[i][j] - sum;
        }
        #pragma parallel for num_threads(number_of_threads)
        for (i = j; i < n; i++) {
            printf("Y\n");
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
}	*/

void crout1(double **A, double **L, double **U, int n, int number_of_threads) {
    int i, j, k;
    double sum = 0;
    #pragma parallel for num_threads(number_of_threads)
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (j = 0; j < n; j++) {
        sum =0 ;
        for (k = 0; k < j; k++) {
            sum = sum + L[j][k] * U[k][j];	
        }
        L[j][j] = A[j][j] - sum;

        #pragma parallel for num_threads(number_of_threads)
        for (i = j; i < n; i++) {
            if(i>j){
                sum = 0;
                for (k = 0; k < j; k++) {
                    sum = sum + L[i][k] * U[k][j];	
                }
                L[i][j] = A[i][j] - sum;
            }
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

void crout2(double **A, double **L, double **U, int n, int number_of_threads) {
    int i,j,k;
    double sum = 0, sum1 = 0, sum2 = 0;
    #pragma omp parallel sections num_threads(number_of_threads)
    {
        #pragma omp section
        {
            for (i = 0; i < n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = n/16; i < 2*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 2*n/16; i < 3*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 3*n/16; i < 4*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 4*n/16; i < 5*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 5*n/16; i < 6*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 6*n/16; i < 7*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 7*n/16; i < 8*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 8*n/16; i < 9*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 9*n/16; i < 10*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 10*n/16; i < 11*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 11*n/16; i < 12*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 12*n/16; i < 13*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 13*n/16; i < 14*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 14*n/16; i < 15*n/16; i++) U[i][i] = 1;
        }

        #pragma omp section
        {
            for (i = 15*n/16; i < n; i++) U[i][i] = 1;
        }
    }
    
    for (j = 0; j < n; j++) {
        sum = 0 ;
        for (k = 0; k < j; k++) {
            sum = sum + L[j][k] * U[k][j];	
        }
        L[j][j] = A[j][j] - sum;
        
        #pragma omp parallel sections num_threads(number_of_threads) private(sum,i,k)
        {
            #pragma omp section
            {
                // printf("A");
                for (i = j+1 + (0*(n-j-1))/8; i < j+1 + (1*(n-j-1))/8; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                // printf("B");
                for (i = j+1 + (1*(n-j-1))/8; i < j+1 + (2*(n-j-1))/8; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                // printf("C");
                for (i = j+1 + (2*(n-j-1)/8); i < j+1 + (3*(n-j-1))/8; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                //printf("D");
                for (i = j+1 + (3*(n-j-1))/8; i < j+1 + (4*(n-j-1))/8; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                //printf("E");
                for (i = j+1 + (4*(n-j-1))/8; i < j+1 + (5*(n-j-1))/8; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                //printf("F");
                for (i = j+1 + (5*(n-j-1))/8; i < j+1 + (6*(n-j-1))/8; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                //printf("G");
                for (i = j+1 + (6*(n-j-1))/8; i < j+1 + (7*(n-j-1))/8; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                //printf("H");
                for (i = j+1 + (7*(n-j-1))/8; i < j+1 + (8*(n-j-1))/8; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                //printf("I");
                for (i = j + (0*(n-j))/8; i < j + (1*(n-j))/8; i++) {
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
            

            #pragma omp section
            {
                //printf("J");
                for (i = j + (1*(n-j))/8; i < j + (2*(n-j))/8; i++) {
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

            #pragma omp section
            {
                //printf("K");
                for (i = j + (2*(n-j))/8; i < j + (3*(n-j))/8; i++) {
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

            #pragma omp section
            {
                //printf("L");
                for (i = j + (3*(n-j))/8; i < j + (4*(n-j))/8; i++) {
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

            #pragma omp section
            {
                //printf("M");
                for (i = j + (4*(n-j))/8; i < j + (5*(n-j))/8; i++) {
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

            #pragma omp section
            {
                //printf("N");
                for (i = j + (5*(n-j))/8; i < j + (6*(n-j))/8; i++) {
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

            #pragma omp section
            {
                //printf("O");
                for (i = j + (6*(n-j))/8; i < j + (7*(n-j))/8; i++) {
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

            #pragma omp section
            {
                //printf("P");
                for (i = j + (7*(n-j))/8; i < j + (8*(n-j))/8; i++) {
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
    }
}	

void crout3(double **A, double **L, double **U, int n, int number_of_threads) {
    int i, j, k;
	double sum = 0;
    #pragma parallel for num_threads(number_of_threads)
	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	for (j = 0; j < n; j++) {
        sum = 0;
        for (k = 0; k < j; k++) {
            sum = sum + L[j][k] * U[k][j];	
        }
        L[j][j] = A[j][j] - sum;
        #pragma omp parallel sections private(i,k,sum)
        {
            #pragma omp section
            {
                #pragma parallel for
                for (i = j; i < n; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];	
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }

            #pragma omp section
            {
                #pragma parallel for
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
	}
}

void print_array(double **A, int n,char *c){ 
    printf("\n--------------------------- %s ---------------------------\n",c);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++)
        printf("%f\t",A[i][j]);
        printf("\n");
    }
    printf("----------------------------------------------------------\n\n");
}

void write_output(char fname[], double** arr, int n){
    FILE *f = fopen(fname, "w");
    for( int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            fprintf(f, "%0.12f ", arr[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

int main(int argc, char* argv[])
{
    int n;
    char *input_filename;
    int strategy;
    double **a;
   
    if (argc != 5) {printf("Enter 5 arguments only. Eg.\"executable_filename n input_filename number_of_threads strategy\"\n");return 0;}
   
    n = atoi(argv[1]);
    input_filename = argv[2];
    number_of_threads = atoi(argv[3]);
    strategy = atoi(argv[4]);
    char *not = argv[3];
   
    a=malloc(n*sizeof(double*)); 
    for(int i=0;i<n;++i)
        a[i]=malloc(n*sizeof(double));
   
    FILE *fptr;
    if ((fptr = fopen(input_filename, "r")) == NULL) {
        printf("Error! can't open the input file.");
        exit(1);
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
        if (!fscanf(fptr, "%lf", &a[i][j])) 
            break;
        }
    }
    fclose(fptr);
   //printf("%d %s %d %d\n",n,input_filename,number_of_threads,strategy);
   
    double **L;
    double **U;
   
    L=malloc(n*sizeof(double*));
    U=malloc(n*sizeof(double*)); 
    for(int i=0;i<n;++i) {
        L[i]=malloc(n*sizeof(double));
        U[i]=malloc(n*sizeof(double));
    }
    if (strategy==0) {
        // printf("X");
        crout0(a,L,U,n);
        // printf("Y");
        char *lf = "output_L_0_";
        char *uf = "output_U_0_";
        char *extension = ".txt";
        char *f1 = (char *) malloc(256);
        char *f2 = (char *) malloc(256);
        strcat(f1,lf); strcat(f1,not); strcat(f1,extension);
        strcat(f2,uf); strcat(f2,not); strcat(f2,extension);

        // printf("|%s|%s|\n",f1,f2);
        write_output(f1,L,n);
        write_output(f2,U,n);
        
        // print_array(L,n,"L");
        // print_array(U,n,"U");
    }
    else if(strategy==1){
        crout1(a,L,U,n,number_of_threads);
        char *lf = "output_L_1_";
        char *uf = "output_U_1_";
        char *extension = ".txt";
        char *f1 = (char *) malloc(256);
        char *f2 = (char *) malloc(256);
        strcat(f1,lf); strcat(f1,not); strcat(f1,extension);
        strcat(f2,uf); strcat(f2,not); strcat(f2,extension);

        write_output(f1,L,n);
        write_output(f2,U,n);

        // print_array(L,n,"L");
        // print_array(U,n,"U");
    }
    else if(strategy==2){
        crout2(a,L,U,n,number_of_threads);
        char *lf = "output_L_2_";
        char *uf = "output_U_2_";
        char *extension = ".txt";
        char *f1 = (char *) malloc(256);
        char *f2 = (char *) malloc(256);
        strcat(f1,lf); strcat(f1,not); strcat(f1,extension);
        strcat(f2,uf); strcat(f2,not); strcat(f2,extension);

        write_output(f1,L,n);
        write_output(f2,U,n);

        // print_array(L,n,"L");
        // print_array(U,n,"U");
    }
    else if(strategy==3){
        crout3(a,L,U,n,number_of_threads);
        char *lf = "output_L_3_";
        char *uf = "output_U_3_";
        char *extension = ".txt";
        char *f1 = (char *) malloc(256);
        char *f2 = (char *) malloc(256);
        strcat(f1,lf); strcat(f1,not); strcat(f1,extension);
        strcat(f2,uf); strcat(f2,not); strcat(f2,extension);

        write_output(f1,L,n);
        write_output(f2,U,n);

        // print_array(L,n,"L");
        // print_array(U,n,"U");
    }
}