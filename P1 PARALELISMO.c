#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define MASTER 0

int main(int argc, char *argv[]) {
    int rank, numprocs, i, done = 0, n;
    double PI25DT = 3.141592653589793238462643;
    double pi, h, sum, x, local_pi;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    while (!done) {
        if (rank == MASTER) {
            printf("Enter the number of intervals: (0 quits)\n");
            scanf("%d", &n);
            for(int j = 1; j<numprocs; j++){
            MPI_Send(&n, 1, MPI_INT, j, 0, MPI_COMM_WORLD);            } 
        }else{
        MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        if (n == 0) break;

        h = 1.0 / (double) n;
        sum = 0.0;
        for (i = rank + 1; i <= n; i += numprocs) {
            x = h * ((double)i - 0.5);
            sum += 4.0 / (1.0 + x*x);
        }
        local_pi = h * sum;

        if (rank != MASTER) {
            MPI_Send(&local_pi, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
        } else {
            pi = local_pi;
            for (i = 1; i < numprocs; i++) {
                MPI_Recv(&local_pi, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                pi += local_pi;
            }
            printf("pi is approx. %.16f, Error: %.16f\n", pi, fabs(pi - PI25DT));
        }
    }

    MPI_Finalize();
    return 0;
}
