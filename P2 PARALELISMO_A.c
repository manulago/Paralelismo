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
        }
        
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (n == 0) break;

        h = 1.0 / (double) n;
        sum = 0.0;
        for (i = rank + 1; i <= n; i += numprocs) {
            x = h * ((double)i - 0.5);
            sum += 4.0 / (1.0 + x*x);
        }
        local_pi = h * sum;

        MPI_Reduce(&local_pi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == MASTER) {
            printf("pi is approx. %.16f, Error: %.16f\n", pi, fabs(pi - PI25DT));
        }
    }

    MPI_Finalize();
    return 0;
}
