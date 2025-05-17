/* cs.cpp */

//
// Performs the actual contrast stretching.
//
// << YOUR NAME >>
//
// Initial author:
//   Prof. Joe Hummel
//   Northwestern University
//

#include "app.h"
#include "matrix.h"
#include <mpi.h>


//
// Given 9 values (imagine a pixel and its 8 surrounding neighbors), returns the median
// value while also returning min and max values (via last 2 parameters)
//
uchar median(uchar a, uchar b, uchar c, uchar d, uchar e, 
						 uchar f, uchar g, uchar h, uchar i, uchar &min, uchar &max)
{
	uchar A[9];

	A[0] = a;  // capture 9 values:
	A[1] = b;
	A[2] = c;
	A[3] = d;
	A[4] = e;
	A[5] = f;
	A[6] = g;
	A[7] = h;
	A[8] = i;

	// selection sort:
	for(int i=0; i < 8; i++)
	{
		int min = i;

		// select smallest in subset:
		for (int j=i+1; j < 9; j++)
		{
			if (A[j] < A[min])
				min = j;
		}

		// swap:
		uchar temp = A[i];
		A[i] = A[min];
		A[min] = temp;
	}//for

	// return min, max & median of 9 elements:
	min = A[0];
	max = A[8];

	return A[4];
}


//
// Given a pixel value (a single RGB value) P and it's 8 neighbors, returns the new
// pixel value P' based on min & max values of neighbors.
//
uchar NewPixelValue(uchar UL, uchar UP, uchar UR,
										uchar L,  uchar P,  uchar R,
										uchar DL, uchar DW, uchar DR, 
	                  int stepby)
{
	uchar newp, med, min, max;
	double ratio;

	med = median(UL, UP, UR, L, P, R, DL, DW, DR, min, max);

	if (min == max)  // pixels are all the same:
		newp = min;
	else
	{
		// are neighbors darker than pixel, or lighter?
		ratio = (P - min) / (double) (max - min);

		if (ratio < 0.5)  // darker, so make pixel darker if possible:
		{
			if (P > stepby) // beware of underflow:
				newp = P - stepby;
			else  
				newp = 0;
		}
		else if (ratio > 0.5)  // lighter, so make pixel lighter if possible:
		{
			if (P < 255-stepby)  // beware of overflow:
				newp = P + stepby;
			else  
				newp = 255;
		}
		else  // neither, so leave pixel alone:
			newp = P;
	}

	// done:
	return newp;
}


//
// Copy over boundary rows to 2nd image matrix:
//
void copy_boundary(uchar** image2, uchar** image, int rows, int cols)
{
	int basecol = 0;

	for (int col = 0; col < cols; col++, basecol += 3)
	{
		image2[0][basecol] = image[0][basecol];
		image2[0][basecol + 1] = image[0][basecol + 1];
		image2[0][basecol + 2] = image[0][basecol + 2];
	}

	basecol = 0;

	for (int col = 0; col < cols; col++, basecol += 3)
	{
		image2[rows - 1][basecol] = image[rows - 1][basecol];
		image2[rows - 1][basecol + 1] = image[rows - 1][basecol + 1];
		image2[rows - 1][basecol + 2] = image[rows - 1][basecol + 2];
	}

	for (int row = 0; row < rows; row++)
	{
		image2[row][0] = image[row][0];
		image2[row][1] = image[row][1];
		image2[row][2] = image[row][2];
	}

	for (int row = 0; row < rows; row++)
	{
		image2[row][cols * 3 - 3] = image[row][cols * 3 - 3];
		image2[row][cols * 3 - 2] = image[row][cols * 3 - 2];
		image2[row][cols * 3 - 1] = image[row][cols * 3 - 1];
	}

}


//
// stretch_one_pixel:
//
void stretch_one_pixel(uchar** image2, uchar** image, int baserow, int basecol)
{
	int prevrow = baserow - 1;  // row above
	int nextrow = baserow + 1;  // row below

	// columns are a little trickier, since a "column" is really 3 physical cols: RGB
	int prevcol = basecol - 3;  // previous column start:
	int nextcol = basecol + 3;  // next column start:

	//
	// now update each pixel based on nearest neighbors around that pixel:
	//

	// Blue:
	image2[baserow][basecol] = NewPixelValue(image[prevrow][prevcol],
		image[prevrow][basecol],
		image[prevrow][nextcol],
		image[baserow][prevcol],
		image[baserow][basecol],
		image[baserow][nextcol],
		image[nextrow][prevcol],
		image[nextrow][basecol],
		image[nextrow][nextcol],
		1 /*stepby*/);

	// Green:
	basecol++;
	prevcol++;
	nextcol++;

	image2[baserow][basecol] = NewPixelValue(image[prevrow][prevcol],
		image[prevrow][basecol],
		image[prevrow][nextcol],
		image[baserow][prevcol],
		image[baserow][basecol],
		image[baserow][nextcol],
		image[nextrow][prevcol],
		image[nextrow][basecol],
		image[nextrow][nextcol],
		1 /*stepby*/);

	// Red:
	basecol++;
	prevcol++;
	nextcol++;

	image2[baserow][basecol] = NewPixelValue(image[prevrow][prevcol],
		image[prevrow][basecol],
		image[prevrow][nextcol],
		image[baserow][prevcol],
		image[baserow][basecol],
		image[baserow][nextcol],
		image[nextrow][prevcol],
		image[nextrow][basecol],
		image[nextrow][nextcol],
		1 /*stepby*/);
}


//
// Performs contrast stetch over the given image --- i.e. makes lighter colors lighter and
// darker colors darker.
//
// NOTE: the image matrix is a 2D matrix of pixels.  However, each pixel is 3 colors RGB,
// which means a "column" in the matrix is actually 3 columns: Blue, Green, and Red.  So
// as we loop over the matrix, we loop by 3 as we go along the column.
//
uchar **ContrastStretch(uchar **image, int rows, int cols, int steps)
{
	//
	// First, we need a temporary 2D matrix of the same size:
	//
	uchar **image2 = New2dMatrix<uchar>(rows, cols*3);

	//
	// copy the boundary rows into the temp image so we can swap
	// pointers after each step:
	//
	copy_boundary(image2, image, rows, cols);

	//
	// Okay, now perform contrast stretching, one step at a time:
	//
	int step = 1;

	while (step <= steps)
	{
		// cout << "** Step " << step << "..." << endl;

		//
		// Okay, for each row (except boundary rows), lighten/darken pixel:
		//
		for (int row = 1; row < rows-1; row++)
		{
			//
			// And for each column (except boundary columns), lighten/darken pixel:
			//
			// columns are a little trickier, since a "column" is really 
			// 3 physical cols: RGB
			//
			int basecol = 3;  // start of column (skip boundary column 0):

			for (int col = 1; col < cols-1; col++, basecol += 3)
			{
				stretch_one_pixel(image2, image, row, basecol);
			}
		}

		//
		// flip the image pointers and step:
		//
		uchar** tempi = image;
		image = image2;
		image2 = tempi;

		step++;

	}//while-each-step

	//
	// done!
	//
	Delete2dMatrix(image2);

	return image;
}

uchar** main_process(uchar** image, int rows, int cols, int steps, int numProcs) {
	cout << "main starting..." << endl;
  
	// approximate number of rows if evenly split among workers (not including boundary rows)
	int rows_per_process = rows / numProcs;
	
	// Send contrast sketch parameters to workers:
	for (int w = 1; w < numProcs; w++) {
		int dest = w;
		int count = 1;
		int tag = 0;
		
		// overlapping boundary rows are sent as well since contrast sketching needs info on boundaries too
		int startRow = (rows_per_process * w) - 1;
		// cout << "startRow: " << startRow << endl;
		
		int endRow;
		if (w == numProcs - 1) {
			// last worker completes the rest of the rows
			endRow = rows;
		}
		else {
			// if not last row, add 2 boundary rows to the chunk
			endRow = startRow + rows_per_process + 2;
		}
		// cout << "endRow: " << endRow << endl;
		// chunk rows is split size + number of boundary rows, aka number of rows in sent chunk
		int chunk_rows = endRow - startRow;

		// first we send the number of rows and columns of the sent image chunk and also the number of steps
		MPI_Send(&steps, count, MPI_INT, dest, tag, MPI_COMM_WORLD);
		// cout << "main sent steps to worker " << w << endl;
		MPI_Send(&chunk_rows, count, MPI_INT, dest, tag, MPI_COMM_WORLD);
		// cout << "main sent chunk rows to worker " << w << endl;
		MPI_Send(&cols, count, MPI_INT, dest, tag, MPI_COMM_WORLD);		// cols is same as main image since we split by rows only
		// cout << "main sent cols to worker " << w << endl;
	
		// now we send the image
		count = chunk_rows * (cols * 3);
		MPI_Send(image[startRow], count, MPI_UNSIGNED_CHAR, dest, tag, MPI_COMM_WORLD);
		// cout << "main sent image chunk to worker " << w << endl;
	}
  
	uchar** result_image = New2dMatrix<uchar>(rows, cols * 3);

	// instead of doing nothing, the main process also performs it's own share of contrast sketching	
	int chunk_rows = min(rows_per_process + 1, rows);	// main process has 1 boundary row, makes sure we don't go out of bounds
	uchar** image_chunk = ContrastStretch(image, chunk_rows, cols, 1);	// will overwrite the image

	// copy the chunk of results to the result image
	for (int i = 0; i < rows_per_process; i++) {
		for (int j = 0; j < cols * 3; j++) {
			result_image[i][j] = image_chunk[i][j];
		}
	}
  
	// now receive the remaining results from the workers:
	cout << "** Receiving results..." << endl;
	cout.flush();
	for (int w = 1; w < numProcs; w++) {
		int startRow = rows_per_process * w;
		int endRow = startRow + rows_per_process;
		
		// last worker completes the rest of the rows
		if (w == numProcs - 1) {
			endRow = rows;
			rows_per_process = endRow - startRow;
		}

		// receive chunk of results from worker:
		int src = w;
		int count = rows_per_process * (cols * 3);
		int tag = 0;

		MPI_Recv(result_image[startRow], count, MPI_UNSIGNED_CHAR, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	// cout << "main process received results from workers" << endl;

	// return pointer to result image
	return result_image;
}

void worker_process(int myRank, int numProcs) {
	cout << "worker " << myRank << " starting..." << endl;
	cout.flush();

	int src = 0;   // main process rank
	int count = 1;  
	int tag = 0;
	int steps, chunk_rows, cols;	// parameters to be received from main
	int i = 0;
	
	while (i < steps) {
		MPI_Recv(&steps, count, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// cout << "worker " << myRank << " received steps from main" << endl;
		MPI_Recv(&chunk_rows, count, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// cout << "worker " << myRank << " received chunk rows from main" << endl;
		MPI_Recv(&cols, count, MPI_INT, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// cout << "worker " << myRank << " received cols from main" << endl;

		// receive the image chunk from main with boundary rows
		uchar** image_chunk = New2dMatrix<uchar>(chunk_rows, cols * 3);
		count = chunk_rows * (cols * 3);	// number of rows in chunk sent from main
		MPI_Recv(image_chunk[0], count, MPI_UNSIGNED_CHAR, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// cout << "worker " << myRank << " received image chunk from main" << endl;
	
		// perform contrast stretching on the image chunk
		uchar** stretched_image = ContrastStretch(image_chunk, chunk_rows, cols, 1);
	
		// send the result image chunk back to main
		int dest = 0;
	
		// calculates size of chunk sent back to main
		if (myRank == numProcs - 1) {
			chunk_rows = chunk_rows - 1; // last worker only has one boundary row
			// cout << "last worker sent chunk rows: " << chunk_rows << endl;
		}
		else {
			chunk_rows = chunk_rows - 2; // all other workers have two boundary rows that we don't want to send back
		}
		count = chunk_rows * (cols * 3);
		
		// skips first row since it is always a boundary row in the worker process (main process handles case where theres no boundary row at start)
		MPI_Send(stretched_image[1], count, MPI_UNSIGNED_CHAR, dest, tag, MPI_COMM_WORLD);
		// cout << "worker " << myRank << " sent results to main for step " << i << endl;
		i++;
	}
}