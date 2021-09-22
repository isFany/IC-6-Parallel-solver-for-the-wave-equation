#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <ctime>
#include <chrono>
#include "Matrix.h"
#include "Solver.h"


int id, p;
int main(int argc, char* argv[]) {

	Grid test_grid(imax - 2, jmax - 2, false, true, false);
	strip test_strip(imax - 2, jmax - 2, false, true, false);
	double t = 0.0;
	double t_out = 0;
	// clock_t start = clock();
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);  
	MPI_Comm_size(MPI_COMM_WORLD, &p); 
	srand(time(NULL) + id * 10);
	test_strip.strip_solver(id, p,t,t_out);
	test_grid.grid_solver(id, p,t,t_out);
	MPI_Finalize();

	// clock_t end = clock();
	// cout << "cost" << (double)(end - start) / CLOCKS_PER_SEC << "seconds" << endl;
}