#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
	int myRank, numProcs;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	if (myRank > 0) { // worker:

		cout << "Hello from worker: rank=" << myRank << endl;

	}
	else { // main:

		cout << "Hello from main: rank=" << myRank
			<< ", total of " << numProcs << " processes"
			<< endl;

	}

	MPI_Finalize();
	return 0;
}
