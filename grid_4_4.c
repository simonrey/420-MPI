//Simon Reynders
//ECSE 420 grid_4_4.c
//Lab 2


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define TAG 0
#define eta 0.0002
#define rho 0.5
#define G 0.75

float interiorCalculate(float * u_array);
int main(int argc, char * argv[]){
	
	// Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //Make a world_group to put all ranks in
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD,&world_group);


    //Make a com between rows
    int color = world_rank / 4;
    MPI_Comm row_comms;
    MPI_Comm_split(MPI_COMM_WORLD,color,world_rank,&row_comms);

    //get the ranks within rows and their size
    int row_rank,row_size;
    MPI_Comm_rank(row_comms,&row_rank);
    MPI_Comm_size(row_comms,&row_size);


    //make a com between columns
    color = (world_rank + 4) % 4;
    MPI_Comm column_comms;
    MPI_Comm_split(MPI_COMM_WORLD,color,world_rank,&column_comms);

    //get the ranks within columns and their size
    int column_rank,column_size;
    MPI_Comm_rank(column_comms,&column_rank);
    MPI_Comm_size(column_comms,&column_size);

    /*
		It helps to imagine the 4x4 arrangement as follows

		|00| - |01| - |02| - |03|
		-------------------------
		|04| - |05| - |06| - |07|
		-------------------------
		|08| - |09| - |10| - |11|
		-------------------------
		|12| - |13| - |14| - |15|

		where 00 has row_rank = 0 and column_rank = 0
		where 01 has row_rank = 1 and column_rank = 0
		where 04 has row_rank = 0 and column_rank = 1
		where 05 has row_rank = 1 and column_rank = 1
		etc

		There are communicators between rows and columns
    */


    //init. all the u's and parse the argument
    float u=0,u1=0,u2=0;
	MPI_Request request;
	int iterations=strtol(argv[1],NULL,10);	
	float u_array[6];

	//adding 1 to u1(N/2,N/2) 
	if(world_rank==10){
		u1+=1;
	}
	FILE * fp;
	if(world_rank==10){
		fp = fopen("prog_output2.h","w");
	}


	MPI_Barrier(MPI_COMM_WORLD);

	for(int i=0;i<iterations;i++){
		int count=1;
		//message passing between inner nodes
		//5,6,9,10,1,2,5,7,8,11,13,14
	    if(row_rank==1 || row_rank==2 || column_rank==1 || column_rank==2){
		    
		    MPI_Status status;
	    	u_array[0]=u1;
	    	u_array[5]=u2;
	    	//for processes with rank 5 and 10 to exchange u1 with internal neighbors and receive u1 from bordering outer ranks
	    	if(world_rank%5==0){
	    		MPI_Sendrecv(&u1,1,MPI_FLOAT,column_rank==1?2:1,TAG,&u_array[count],1,MPI_FLOAT,column_rank==1?2:1,TAG,column_comms,&status);
		   		count++;
		   		MPI_Sendrecv(&u1,1,MPI_FLOAT,row_rank==1?2:1,TAG,&u_array[count],1,MPI_FLOAT,row_rank==1?2:1,TAG,row_comms,&status);
		   		count++;
		   		MPI_Recv(&u_array[count],1,MPI_FLOAT,column_rank==1?0:3,TAG,column_comms,MPI_STATUS_IGNORE);
		   		count++;
				MPI_Recv(&u_array[count],1,MPI_FLOAT,row_rank==1?0:3,TAG,row_comms,MPI_STATUS_IGNORE);
	    	
	    		//once every internal ranks has their u1's from their neighbors, calculate u's
	    		u=interiorCalculate(u_array);
	    	}
	    	//for processes with rank 6 and 9 to exchange u1 with internal neighbors and receive u1 from bordering outer ranks
	    	else if(world_rank%3==0){
	    		MPI_Sendrecv(&u1,1,MPI_FLOAT,column_rank==1?2:1,TAG,&u_array[count],1,MPI_FLOAT,column_rank==1?2:1,TAG,column_comms,&status);
	    		count++;
		   		MPI_Sendrecv(&u1,1,MPI_FLOAT,row_rank==1?2:1,TAG,&u_array[count],1,MPI_FLOAT,row_rank==1?2:1,TAG,row_comms,&status);
	    		count++;
		   		MPI_Recv(&u_array[count],1,MPI_FLOAT,column_rank==1?0:3,TAG,column_comms,MPI_STATUS_IGNORE);
				count++;
				MPI_Recv(&u_array[count],1,MPI_FLOAT,row_rank==1?0:3,TAG,row_comms,MPI_STATUS_IGNORE);
				
				//once every internal ranks has their u1's from their neighbors, calculate u's
	    		u=interiorCalculate(u_array);
	    	}
	    	//for bordering outer ranks - 1,2,13,14 - to send u1 to their inner neighbors
	    	else if(row_rank==1||row_rank==2){
	    		MPI_Send(&u1,1,MPI_FLOAT,column_rank==0?1:2,TAG,column_comms);
	    	}
	    	//for bordering outer ranks - 4,7,8,11 - to send u1 to their inner neighbors
	    	else{
	    		MPI_Send(&u1,1,MPI_FLOAT,row_rank==0?1:2,TAG,row_comms);
	    	}
	    }

	    //Wait here until all the internal ranks are ready to move forward
	    MPI_Barrier(MPI_COMM_WORLD);

	    //captures all the sides
	    //1,2,4,7,8,11,13,14
	    if(((row_rank==1 || row_rank ==2) && (column_rank==0 || column_rank==3)) || ((row_rank==0 || row_rank==3) && (column_rank==1 || column_rank==2))){
	    	// printf("%i\n", world_rank);

	    	//side ranks - 1,2,13,14 - receive u from their internal neighbors
	    	if(world_rank<4 || world_rank>12){
	    		MPI_Recv(&u,1,MPI_FLOAT,column_rank==0?1:2,TAG,column_comms,MPI_STATUS_IGNORE);
		    	u = G * u;
	    	}
	    	//side ranks - 4,7,8,11 - receive u from their internal neighbors
	    	else{
	    		MPI_Recv(&u,1,MPI_FLOAT,row_rank==0?1:2,TAG,row_comms,MPI_STATUS_IGNORE);
	    		u = G * u;
	    	}
	    }
	    //captures internal ranks - 5,6,9,10 - to send their NEW u's as calculated in the interal step to their neighboring side ranks
		if(row_rank==1 || row_rank==2 || column_rank==1 || column_rank==2){
		   	if(world_rank%5==0){
		   		MPI_Send(&u,1,MPI_FLOAT,column_rank==1?0:3,TAG,column_comms);
		   		MPI_Send(&u,1,MPI_FLOAT,row_rank==1?0:3,TAG,row_comms);
		   	}
		   	else if(world_rank%3==0){
				MPI_Send(&u,1,MPI_FLOAT,column_rank==1?0:3,TAG,column_comms);
		   		MPI_Send(&u,1,MPI_FLOAT,row_rank==1?0:3,TAG,row_comms);
		   	}
		}

		//Wait here until all the message passing between the internal ranks and their side neighbors is complete
		MPI_Barrier(MPI_COMM_WORLD);

		//captures all the corners, has those ranks receive from their neighbors underneath
		//0,3,12,15
		if((row_rank==0||column_rank==0||row_rank==3||column_rank==3) && world_rank%3==0){
			MPI_Recv(&u,1,MPI_FLOAT,column_rank==0?1:2,TAG,column_comms,MPI_STATUS_IGNORE);
			// MPI_Recv(&u,1,MPI_FLOAT,row_rank==0?1:2,TAG,row_comms,MPI_STATUS_IGNORE);
			u = G * u;
		}
		//captures all the ones sending to above
		//4,7,8,11
		if((row_rank==0||row_rank==3) && world_rank%3!=0){
			MPI_Send(&u,1,MPI_FLOAT,column_rank==1?0:3,TAG,column_comms);
			// MPI_Send(&u,1,MPI_FLOAT,row_rank==1?0:3,TAG,row_comms);
		}

		//Wait till all ranks get here
		MPI_Barrier(MPI_COMM_WORLD);

		//Update everything for next iteration
		u2=u1;
		u1=u;
		//Only printing (N/2,N/2) which in this code is always world_rank = 10
		if(world_rank==10){
			
            if(i==iterations-1){
                // printf("%f\n", u[N/2][N/2]);
                fprintf(fp,"%f\n",u);
                fclose(fp);
            }
            else{
                // printf("%f,\n", u[N/2][N/2]);
                fprintf(fp,"%f,\n",u);
            }
		}

		//Wait again for the printing
		MPI_Barrier(MPI_COMM_WORLD);
		//GONE THROUGH 1 ITERATION!
	}

	//Wait till everyone gets here before freeing up all the comms that were made
	MPI_Barrier(MPI_COMM_WORLD);

    MPI_Group_free(&world_group);
    MPI_Comm_free(&row_comms);
    MPI_Comm_free(&column_comms);
	MPI_Finalize();

}

float interiorCalculate(float * u_array){
	float t0,t1,t2,t3;
	float numerator, denominator;

	t0 = u_array[1]+u_array[2]+u_array[3]+u_array[4]-(4*u_array[0]);
	t1 = rho*t0;
	t2 = (2*u_array[0]);
	t3 = (1-eta)*u_array[5];

	denominator = 1 + eta;
	numerator = t1 + t2 - t3;

	float result = numerator/denominator;

	return result;
}