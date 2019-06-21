#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// Initialize constants
double lambda;
int N, M;
double omega, tau;

int main(int argc, char* argv[]){
    lambda = 100.0;
    double func_g = 0;

    // Start runtime clock

    clock_t start_time, end_time;
    start_time = clock();

    // Set up args to be taken as inputs

    N = atoi(argv[1]);
    M = 2*N - 1;
    omega = atof(argv[2]);
    tau = atof(argv[3]);

    double delta = 2.0/(N - 1);

    // Initialize matrix

    double res_max = 1.0 + tau;
    double *u[M], *Res[M];
    int i, j;
    for (i = 0; i < M; ++i){
        u[i] = (double *)malloc((N) * sizeof(double));
        Res[i] = (double *)malloc((N) * sizeof(double));
        for (j = 0; j < N; ++j){
            u[i][j] = 0;
            Res[i][j] = 0;
            }
        }

    // Implement SOR algorithm

    // Boundary i = 0
    int cnt = 0;
    while(res_max > tau){
       res_max = 0;
        cnt++;

        // Inside domain now

        for (i = 0; i < M; ++i){
            for(j = 1; j < N-1 ; ++j){
                if(i==0){

                // Forcing term

                func_g = 10*lambda/sqrt(M_PI)*(exp(-lambda*lambda*(9 + (-1 + j*delta)*(-1 + j*delta)))
                        -exp(-lambda*lambda*(1 + (-1 + j*delta)*(-1 + j*delta))));

                // Use ghost point to update matrix

                Res[0][j] = (u[i+1][j] + u[i][j+1] + u[i][j-1]
                            -3*u[i][j])/4 - delta*delta*func_g/4;

                // Add updated step to solution

                u[0][j] = u[0][j] + omega*Res[0][j];

                // Add updated step to residue term

                if(fabs(Res[0][j]) > res_max)
                    res_max = fabs(Res[0][j]);
            }

            else if (i<M-1){

            // Forcing term

			func_g = 10*lambda/sqrt(M_PI)*(exp(-lambda*lambda*((-3 + i*delta)*(-3 + i*delta) + (-1 + j*delta)*(-1 + j*delta)))
                      -exp(-lambda*lambda*((-1 + i*delta)*(-1 + i*delta) + (-1 + j*delta)*(-1 + j*delta))));

            // Use neighboring points to update

			Res[i][j] = -u[i][j] + (u[i][j+1] + u[i][j-1]
                        +u[i-1][j] + u[i+1][j])/4 - delta*delta*func_g/4;

            // Add updated step to solution

			u[i][j] = u[i][j] + omega*Res[i][j];

            // Add updated step to residue

			if(fabs(Res[i][j]) > res_max)
                res_max = fabs(Res[i][j]);
            }

            else if (i=M-1){

        // Forcing term

		func_g = 10*lambda/sqrt(M_PI)*(exp(-lambda*lambda*(1 + (-1 + j*delta)*(-1 + j*delta)))
                        -exp(-lambda*lambda*(9 + (-1 + j*delta)*(-1 + j*delta))));

		//Use ghost point to update matrix

		Res[M-1][j] = (u[i-1][j] + u[i][j+1] + u[i][j-1]
                            -3*u[i][j])/4 - delta*delta*func_g/4;

        // Add updated step to solution

		u[M-1][j] = u[M-1][j] + omega*Res[M-1][j];

        // Add updated step to residue

		if(fabs(Res[M-1][j]) > res_max)
            res_max = fabs(Res[M-1][j]);}

        }
    }

}

// Output iteration to console

printf("Number of iterations: %d\n", cnt);

// Write soln to target file

FILE *fileid_u = fopen("u.out", "w");
for(i = 0; i < M; ++i){
	fwrite(u[i], sizeof(double), N, fileid_u);
}
fclose(fileid_u);

// Free memory

for(i = 0; i < M; ++i){
   free(u[i]);
   free(Res[i]);
}

end_time = clock();

// Output run time to console
printf("Runtime: %1.4fs\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
return 0;
}
