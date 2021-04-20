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


void print_array(double **A, int n,char *c){
    printf("\n--------------------------- %s ---------------------------\n",c);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++)
        printf("%f\t",A[i][j]);
        printf("\n");
    }
    printf("----------------------------------------------------------\n\n");
}

void crout1(double **A, double **L, double **U, int n, int number_of_threads) {
    int i, j, k;
    double sum = 0;
    #pragma parallel for num_threads(number_of_threads)
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }

    print_array(A,n,"A");
    print_array(L,n,"L");
    print_array(U,n,"U");

    #pragma parallel for num_threads(number_of_threads)
    for (j = 0; j < n; j++) {

    
            for (i = j; i < n; i++) {
                printf("X\n");
                sum = 0;
                for (k = 0; k < j; k++) {
                    sum = sum + L[i][k] * U[k][j];	
                }
                L[i][j] = A[i][j] - sum;
                // printf("j:%d\ti:%d\n",j,i);
                // print_array(L,n,"L");
            }
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
        crout0(a,L,U,n);
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
}