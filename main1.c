#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>

int number_of_threads;

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
/*
void crout0(double const **A, double **L, double **U, int n) {
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
*/

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
}	
*/

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

void strategy2(double **A, double **L, double **U, int n, int t) {
	#pragma omp parallel sections num_threads(t)
	{
		#pragma omp section
		{
			for (int i = (n*0)/16; i < (n*1)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*1)/16; i < (n*2)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*2)/16; i < (n*3)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*3)/16; i < (n*4)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*4)/16; i < (n*5)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*5)/16; i < (n*6)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*6)/16; i < (n*7)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*7)/16; i < (n*8)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*8)/16; i < (n*9)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*9)/16; i < (n*10)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*10)/16; i < (n*11)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*11)/16; i < (n*12)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*12)/16; i < (n*13)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*13)/16; i < (n*14)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*14)/16; i < (n*15)/16; i++) {
				U[i][i] = 1;
			}
		}
		#pragma omp section
		{
			for (int i = (n*15)/15; i < (n*16)/16; i++) {
				U[i][i] = 1;
			}
		}
	}
	for (int j = 0; j < n; j++) {
		double sum = 0;
		for (int k = 0; k < j; k++) {
			sum = sum + L[j][k] * U[k][j];	
		}
		L[j][j] = A[j][j] - sum;
		#pragma omp parallel sections num_threads(t)
		{
			#pragma omp section
			{
				for (int i = (j+1) + (0*(n-j-1))/8; i < (j+1) + (1*(n-j-1))/8; i++) {
					int sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = (j+1) + (1*(n-j-1))/8; i < (j+1) + (2*(n-j-1))/8; i++) {
					int sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = (j+1) + (2*(n-j-1))/8; i < (j+1) + (3*(n-j-1))/8; i++) {
					int sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = (j+1) + (3*(n-j-1))/8; i < (j+1) + (4*(n-j-1))/8; i++) {
					int sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = (j+1) + (4*(n-j-1))/8; i < (j+1) + (5*(n-j-1))/8; i++) {
					int sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = (j+1) + (5*(n-j-1))/8; i < (j+1) + (6*(n-j-1))/8; i++) {
					int sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = (j+1) + (6*(n-j-1))/8; i < (j+1) + (7*(n-j-1))/8; i++) {
					int sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = (j+1) + (7*(n-j-1))/8; i < (j+1) + (8*(n-j-1))/8; i++) {
					int sum = 0;
					for (int k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];	
					}
					L[i][j] = A[i][j] - sum;
				}
			}
			#pragma omp section
			{
				for (int i = j + (0*(n-j))/8; i < j + (1*(n-j))/8; i++) {
					int sum = 0;
					for(int k = 0; k < j; k++) {
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
				for (int i = j + (1*(n-j))/8; i < j + (2*(n-j))/8; i++) {
					int sum = 0;
					for(int k = 0; k < j; k++) {
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
				for (int i = j + (2*(n-j))/8; i < j + (3*(n-j))/8; i++) {
					int sum = 0;
					for(int k = 0; k < j; k++) {
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
				for (int i = j + (3*(n-j))/8; i < j + (4*(n-j))/8; i++) {
					int sum = 0;
					for(int k = 0; k < j; k++) {
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
				for (int i = j + (4*(n-j))/8; i < j + (5*(n-j))/8; i++) {
					int sum = 0;
					for(int k = 0; k < j; k++) {
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
				for (int i = j + (5*(n-j))/8; i < j + (6*(n-j))/8; i++) {
					int sum = 0;
					for(int k = 0; k < j; k++) {
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
				for (int i = j + (6*(n-j))/8; i < j + (7*(n-j))/8; i++) {
					int sum = 0;
					for(int k = 0; k < j; k++) {
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
				for (int i = j + (7*(n-j))/8; i < j + (8*(n-j))/8; i++) {
					int sum = 0;
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

void crout2(double **A, double **L, double **U, int n, int number_of_threads) {
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
    
    int i,j,k;
    double sum = 0, sum1 = 0, sum2 = 0;
    for (j = 0; j < n; j++) {
        sum = 0 ;
        for (k = 0; k < j; k++) {
            sum = sum + L[j][k] * U[k][j];	
        }
        L[j][j] = A[j][j] - sum;
        
        #pragma omp parallel sections num_threads(number_of_threads)
        {
            #pragma omp section
            {
                printf("A");
                for (i = j+1; i < n/4; i++) {
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

            #pragma omp section
            {
                printf("B");
                for (i = n/4; i < 2*n/4; i++) {
                    if(i>j){
                        sum2 = 0;
                        for (k = 0; k < j; k++) {
                            sum2 = sum2 + L[i][k] * U[k][j];	
                        }
                        L[i][j] = A[i][j] - sum2;
                    }
                    sum2 = 0;
                    for(k = 0; k < j; k++) {
                        sum2 = sum2 + L[j][k] * U[k][i];
                    }
                    if (L[j][j] == 0) {				
                        exit(0);
                    }
                    U[j][i] = (A[j][i] - sum2) / L[j][j];
                }
            }
            
            #pragma omp section
            {
                printf("C");
                for (i = 2*n/4; i < 3*n/4; i++) {
                    if(i>j){
                        sum2 = 0;
                        for (k = 0; k < j; k++) {
                            sum2 = sum2 + L[i][k] * U[k][j];	
                        }
                        L[i][j] = A[i][j] - sum2;
                    }
                    sum2 = 0;
                    for(k = 0; k < j; k++) {
                        sum2 = sum2 + L[j][k] * U[k][i];
                    }
                    if (L[j][j] == 0) {				
                        exit(0);
                    }
                    U[j][i] = (A[j][i] - sum2) / L[j][j];
                }
            }

            #pragma omp section
            {
                printf("D");
                for (i = 3*n/4; i < n; i++) {
                    if(i>j){
                        sum2 = 0;
                        for (k = 0; k < j; k++) {
                            sum2 = sum2 + L[i][k] * U[k][j];	
                        }
                        L[i][j] = A[i][j] - sum2;
                    }
                    sum2 = 0;
                    for(k = 0; k < j; k++) {
                        sum2 = sum2 + L[j][k] * U[k][i];
                    }
                    if (L[j][j] == 0) {				
                        exit(0);
                    }
                    U[j][i] = (A[j][i] - sum2) / L[j][j];
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
    printf("%d\n\n\n",strategy);
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

        print_array(L,n,"L");
        print_array(U,n,"U");
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

        print_array(L,n,"L");
        print_array(U,n,"U");
    }
    else if(strategy==2){
        strategy2(a,L,U,n,number_of_threads);
        char *lf = "output_L_2_";
        char *uf = "output_U_2_";
        char *extension = ".txt";
        char *f1 = (char *) malloc(256);
        char *f2 = (char *) malloc(256);
        strcat(f1,lf); strcat(f1,not); strcat(f1,extension);
        strcat(f2,uf); strcat(f2,not); strcat(f2,extension);

        write_output(f1,L,n);
        write_output(f2,U,n);

        print_array(L,n,"L");
        print_array(U,n,"U");
    }
}