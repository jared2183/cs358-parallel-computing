/* main.cpp */

//
// Parallelizes a generic "work matrix" where work is randomly
// distributed in an NxN matrix. Naive parallelization works,
// but doesn't scale. A much more dynamic solution is needed.
// 
// Usage:
//   work [-?] [-t NumThreads]
//
// Author:
//   << YOUR NAME >>
//   Northwestern University
// 
// Initial template:
//   Prof. Joe Hummel
//   Northwestern University
//

#include <iostream>
#include <string>
#include <cstring>
#include <chrono>
#include <sys/sysinfo.h>

#include "alloc2D.h"
#include "workmatrix.h"

using namespace std;


//
// Globals:
//
static int _numThreads = 1;  // default to sequential execution
static int cells = 0;

//
// Function prototypes:
//
static void ProcessCmdLineArgs(int argc, char* argv[]);


//
// main:
//
int main(int argc, char *argv[])
{
	cout << "** Work Matrix Application **" << endl;
	cout << endl;

	//
	// Set defaults, process environment & cmd-line args:
	//
	ProcessCmdLineArgs(argc, argv);

	WorkMatrix wm;  // NOTE: wm MUST be created in sequential code.

	cout << "Matrix size:  " << wm.num_rows() << "x" << wm.num_cols() << endl;
	cout << "# of threads: " << _numThreads << endl;
	cout << endl;

	cout << "working";
	cout.flush();

	//
	// Solve each cell in the work matrix. Compute time for speedup
	// calculations.
	//
  auto start = chrono::high_resolution_clock::now();
	// saves number of rows and columns to a variable
	int num_rows = wm.num_rows();
	int num_cols = wm.num_cols();

	#pragma omp parallel for schedule(dynamic) num_threads(_numThreads)
	// flattens the array so that dynamic scheduling is done element-by-element instead of row-by-row
	for (int i = 0; i < num_rows * num_cols; i++) {
		// does work in cell [r][c] where r = i / num_rows (integer division) and c = i % num_rows (modulo)
		wm.do_work(i / num_rows, i % num_rows);

		// show some output every 100 cells so we see progress:
		cells++;

		if (cells % 100 == 0) {
			cout << ".";
			cout.flush();
		}
	}
  
  auto stop = chrono::high_resolution_clock::now();
  auto diff = stop - start;
  auto duration = chrono::duration_cast<chrono::milliseconds>(diff);

	cout << endl;
	cout << endl;

  cout << endl;
  cout << "** Done!  Time: " << duration.count() / 1000.0 << " secs" << endl;
	cout << "** Execution complete **" << endl;
  cout << endl;

	return 0;
}


//
// processCmdLineArgs:
//
static void ProcessCmdLineArgs(int argc, char* argv[])
{
	for (int i = 1; i < argc; i++)
	{

		if (strcmp(argv[i], "-?") == 0)  // help:
		{
			cout << "**Usage: work [-?] [-t NumThreads]" << endl << endl;
			exit(0);
		}
		else if ((strcmp(argv[i], "-t") == 0) && (i+1 < argc))  // # of threads:
		{
			i++;
			_numThreads = atoi(argv[i]);
		}
		else  // error: unknown arg
		{
			cout << "**Unknown argument: '" << argv[i] << "'" << endl;
			cout << "**Usage: work [-?] [-t NumThreads]" << endl << endl;
			exit(0);
		}

	}//for
}
