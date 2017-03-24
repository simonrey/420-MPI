//Simon Reynders
//ECSE 420 grid_512_512.c
//Lab 2

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

#define TAG 0
#define eta 0.0002
#define rho 0.5
#define G 0.75
#define N 512


int main(int argc, char * argv[]) {
    /* Number of iterations */
    int iterations = strtol(argv[1], NULL, 0);

    /* Initialize the MPI environment */
    MPI_Init(&argc, &argv);
    /* Get the number of processes */
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //make a_comm
    int color = rank/2;
    MPI_Comm a_comm;
    MPI_Comm_split(MPI_COMM_WORLD,color,rank,&a_comm);
    
    //make b_comm
    color = (rank-1)/2;
    MPI_Comm b_comm;
    MPI_Comm_split(MPI_COMM_WORLD,color,rank,&b_comm);
    

    int a_comm_rank, a_comm_size;
    MPI_Comm_rank(a_comm,&a_comm_rank);
    MPI_Comm_size(a_comm,&a_comm_size);

    int b_comm_rank, b_comm_size;
    MPI_Comm_rank(b_comm,&b_comm_rank);
    MPI_Comm_size(b_comm,&b_comm_size);

    /* Initialize our matrices */
    float u[N][N] = {0};
    float u1[N][N] = {0};
    float u2[N][N] = {0};

    /* Initial condition */
    u1[N/2][N/2] += 1;


    int block_size = N * N / size;
    int num_rows_each = N / size;

    int start = rank * num_rows_each;
    int end = num_rows_each + start -1;

    MPI_Status status;
    FILE * fp;

    if(rank==(size/2)){
        fp = fopen("prog_output.h","w");
    }
    clock_t tic = clock();
    for (int i = 0; i < iterations; i++){

        if(rank==(size-1)){
            MPI_Sendrecv(&u1[start][0],N,MPI_FLOAT,a_comm_rank-1,TAG,&u1[start-1][0],N,MPI_FLOAT,a_comm_rank-1,TAG,a_comm,&status);
        }
        else if(rank==0){
            MPI_Sendrecv(&u1[end][0],N,MPI_FLOAT,a_comm_rank+1,TAG,&u1[end+1][0],N,MPI_FLOAT,a_comm_rank+1,TAG,a_comm,&status);
        }
        else if(rank%2==0){
            MPI_Sendrecv(&u1[end][0],N,MPI_FLOAT,a_comm_rank+1,TAG,&u1[end+1][0],N,MPI_FLOAT,a_comm_rank+1,TAG,a_comm,&status);
            MPI_Sendrecv(&u1[start][0],N,MPI_FLOAT,b_comm_rank-1,TAG,&u1[start-1][0],N,MPI_FLOAT,b_comm_rank-1,TAG,b_comm,&status);
        }
        else{
            MPI_Sendrecv(&u1[start][0],N,MPI_FLOAT,a_comm_rank-1,TAG,&u1[start-1][0],N,MPI_FLOAT,a_comm_rank-1,TAG,a_comm,&status);
            MPI_Sendrecv(&u1[end][0],N,MPI_FLOAT,b_comm_rank+1,TAG,&u1[end+1][0],N,MPI_FLOAT,b_comm_rank+1,TAG,b_comm,&status);   
        }

        // MPI_Barrier(MPI_COMM_WORLD);

        
        /* First we update the center, as shown in the lab instructions */
        for(int i = start; i <= end; i++){
            for(int j = 1; j <= N-2; j++){
                u[i][j] = ((rho) * (u1[i-1][j] + u1[i+1][j] + u1[i][j-1] + u1[i][j+1] - (4 * u1[i][j])) + 2 * u1[i][j] - (1 - eta) * u2[i][j])/(1 + eta);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /* We then update the side elements, which depend on the center elements */
        for(int i = 1; i <= N - 2; i++){

                //sides
                u[0][i] = G * u[1][i];
                u[N-1][i] = G * u[N-2][i];
                u[i][0] = G * u[i][1];
                u[i][N-1] = G * u[i][N-2];
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /* Finally, we now can update the corner elements, which depend on its neighbouring elements */
        u[0][0] = G * u[1][0];
        u[N-1][0] = G * u[N-2][0];
        u[0][N-1] = G * u[0][N-2];
        u[N-1][N-1] = G * u[N-1][N-2];


        // MPI_Barrier(MPI_COMM_WORLD);

        if(rank==(size/2)){
            if(i==iterations-1){
                // printf("%f\n", u[N/2][N/2]);
                fprintf(fp,"%f\n",u[N/2][N/2]);
                fclose(fp);
            }
            else{
                // printf("%f,\n", u[N/2][N/2]);
                fprintf(fp,"%f,\n",u[N/2][N/2]);
            }
        }
        

        /* Here we update, basically copying the arrays*/
        memcpy(u2, u1, sizeof(u2));
        memcpy(u1, u, sizeof(u1));

        MPI_Barrier(MPI_COMM_WORLD);
  
    }
  	clock_t toc = clock();
    
    if(rank == (size/2)){
        fp=fopen("prog_output.h","a");
        fprintf(fp,"Runtime of process with u(N/2,N/2): %f\n",(double)((toc-tic)/CLOCKS_PER_SEC));
        // printf("%f\n",(double)((toc-tic)/CLOCKS_PER_SEC))
    }

    MPI_Comm_free(&a_comm);
    MPI_Comm_free(&b_comm);
    MPI_Finalize();
    return 0;

}