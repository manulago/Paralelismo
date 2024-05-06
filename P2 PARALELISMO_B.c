#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define MASTER 0


int MPI_FlattreeColectiva(void *buff, void *recvbuff, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
    int rank, size;
    double aux, n;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    if(rank != root){
        MPI_Send(buff, count, datatype, root, 0, comm);
    }else{
        aux = *(double*) buff;
        for(int i = 0; i<size; i++){
            if(i != root){
                MPI_Recv(&n, 1, datatype, MPI_ANY_SOURCE, 0, comm, &status);
                aux += n;
            }
        }
        *(double*)recvbuff = aux;
    }
    return status.MPI_ERROR;
}


int MPI_BinomialBcast(void *buff, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
    int numprocesos, rank, proceso,source;

    MPI_Comm_size(comm, &numprocesos);
    MPI_Comm_rank(comm, &rank);

    for (int i = 1; i < numprocesos; i++){
        if(rank < pow(2,i-1)){
            proceso = rank + pow(2,i-1);
            if(proceso < numprocesos){
                printf("Soy el proceso %d y envío a %d\n", rank, proceso);
                MPI_Send(buff, count, datatype, proceso, 0, comm);
            }
        }
        else{
            source =  rank - pow(2,i-1);
            if(source < pow(2,i-1)){
            MPI_Recv(buff,count, datatype,source,0,comm,MPI_STATUS_IGNORE);
            }
        }
    }
    return MPI_SUCCESS;
}



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

        MPI_BinomialBcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (n == 0) break;

        h = 1.0 / (double) n;
        sum = 0.0;
        for (i = rank + 1; i <= n; i += numprocs) {
            x = h * ((double)i - 0.5);
            sum += 4.0 / (1.0 + x*x);
        }
        local_pi = h * sum;

        MPI_FlattreeColectiva(&local_pi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == MASTER) {
            printf("pi is approx. %.16f, Error: %.16f\n", pi, fabs(pi - PI25DT));
        }
    }

    MPI_Finalize();
    exit (0);
}

