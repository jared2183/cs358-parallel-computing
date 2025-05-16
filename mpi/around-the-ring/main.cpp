#include <iostream>
#include <mpi.h>

using namespace std;

void main_process(int numProcs)
{
	MPI_Status status;
	int value;

	cout << "Starting with a total of " << numProcs << " processes" << endl;
	cout << endl;

	cout << "Please enter starting value> ";
	cin >> value;
	cout << endl;

	int dest = 1;  // start by sending to worker 1:
	int count = 1; // 1 value to send
	int tag = 0;   // doesn't matter as long as we all use 0

	MPI_Send(&value, count, MPI_INT, dest, tag, MPI_COMM_WORLD);

	int src = numProcs - 1;  // wait for msg from last worker, can also set to receive from any machine (MPI_ANY_SOURCE)
	MPI_Recv(&value, count, MPI_INT, src, tag, MPI_COMM_WORLD, &status);

	cout << "Main received " << value << endl;

	value++;  // we incr as well

	cout << endl;
	cout << "Main: final value: " << value << endl;
}

void worker_process(int myRank, int numProcs)
{
	MPI_Status status;

	int src = myRank - 1;  // receive from worker to our "left"
	int count = 1;         // 1 value to send
	int tag = 0;           // doesn't matter as long as we all use 0

	int value;

	MPI_Recv(&value, count, MPI_INT, src, tag, MPI_COMM_WORLD, &status);

	cout << "Worker " << myRank << " received " << value << endl;
	cout.flush();

	value++;  // every worker increments

	int dest = myRank + 1;     // send to the worker on our "right"

	if (dest == numProcs)  // if dest doesn't exist, send to main
		dest = 0;

	MPI_Send(&value, count, MPI_INT, dest, tag, MPI_COMM_WORLD);
}

int main(int argc, char* argv[])
{
	int myRank, numProcs;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	if (myRank > 0)
		worker_process(myRank, numProcs);
	else
		main_process(numProcs);

	MPI_Finalize();
	return 0;
}
