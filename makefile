CC=mpicc

all: program1 program2

program1: grid_4_4.c
	$(CC) -o grid_4_4 grid_4_4.c
program2: grid_512_512.c
	$(CC) -o grid_512_512 grid_512_512.c
