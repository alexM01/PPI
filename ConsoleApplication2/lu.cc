#include "stdafx.h"
#include "matrix.h"
using namespace std;


/* Code taken from the GLIBC manual.
*
* Subtract the ‘struct timespec’ values X and Y,
* storing the result in RESULT.
* Return 1 if the difference is negative, otherwise 0.*/

static int
timespec_subtract(struct timespec *result,
	struct timespec *x,
	struct timespec *y)
{
	/* Perform the carry for the later subtraction by updating y.*/
	if (x->tv_nsec < y->tv_nsec) {
		int nsec = (y->tv_nsec - x->tv_nsec) / 1000000000 + 1;
		y->tv_nsec -= 1000000000 * nsec;
		y->tv_sec += nsec;
	}
	if (x->tv_nsec - y->tv_nsec > 1000000000) {
		int nsec = (x->tv_nsec - y->tv_nsec) / 1000000000;
		y->tv_nsec += 1000000000 * nsec;
		y->tv_sec -= nsec;
	}

	/* Compute the time remaining to wait.
	tv_nsec is certainly positive.*/
	result->tv_sec = x->tv_sec - y->tv_sec;
	result->tv_nsec = x->tv_nsec - y->tv_nsec;

	/* Return 1 if result is negative.*/
	return x->tv_sec < y->tv_sec;
}


/* Global variables holding the matrix data. To complete this assignment
* you are requested to only use arrays and access these arrays with
* subscripts. Do not use pointers.*/



#ifdef WIN32
#define CLOCK_REALTIME 0
static int clock_gettime(int, struct timespec *spec)      //C-file part
{
	__int64 wintime; GetSystemTimeAsFileTime((FILETIME*)&wintime);
	wintime -= 116444736000000000i64;  //1jan1601 to 1jan1970
	spec->tv_sec = wintime / 10000000i64;           //seconds
	spec->tv_nsec = wintime % 10000000i64 * 100;      //nano-seconds
	return 0;
}
#endif
const int max_n_elements = 131072;
const int max_n_rows = 16384;

static int n_colsA, n_rowsA;
static double values[max_n_elements];

static int col_ind[max_n_elements];
static int col_indBackup[max_n_elements];
static int row_ptr_begin[max_n_rows];
static int row_ptr_end[max_n_rows];
static int pBsize;
static double valuesL[max_n_elements];
static double valuesU[max_n_elements];
static double valuesLU[max_n_elements];
static double col_indL[max_n_elements];
static double col_indLU[max_n_elements];
static double col_indU[max_n_elements];
static int row_ptr_beginLU[max_n_rows];
static int row_ptr_beginL[max_n_rows];
static int row_ptr_beginU[max_n_rows];
static int row_ptr_beginUBackup[max_n_rows];
static int  row_ptr_endL[max_n_rows];
static int  row_ptr_endLU[max_n_rows];
static int row_ptr_endU[max_n_rows];
static int skip[max_n_rows];
static int skip2[max_n_rows];
static double results[max_n_elements];
static double resultsx[max_n_elements];
static double resultsb[max_n_elements];
static double xvector[max_n_elements];
static double xvector2[max_n_elements];

static int PivotIsZero() {
	//shuffle row of input matrix
	int upper, lower, temp1, temp2;
	bool found = false;
	for (int i = 0; i < pBsize; i++) {
		//upper = col_ind[row_ptr_begin[i + 1]];
		lower = col_ind[row_ptr_begin[i] + i];
		if (lower == i) {
			if (abs(values[row_ptr_begin[i] + i]) < 1) {
				for (int j = i - 1; j >= 0; j++) {
					if (abs(values[row_ptr_begin[j] + j]) > 1) {
						if (!found) {
							found = true;
							temp1 = i; temp2 = j;
						}
						if (row_ptr_end[i] - row_ptr_begin[i] == row_ptr_end[j] - row_ptr_begin[j]) {
							printf("Found");
							return i, j;
						}

					}
				}
			}
		}
	}
	if (found) {
		return temp1, temp2;
	}
	else {
		return -1, -1;
	}
}

static bool switchRows(int currentRow, int newRow) {
	int cBegin, cEnd, nBegin, nEnd, temp;
	if (currentRow > 0) {
		//switch rows currentRow and newRow
		cBegin = row_ptr_begin[currentRow];
		cEnd = row_ptr_end[currentRow];
		nBegin = row_ptr_begin[newRow];
		nEnd = row_ptr_end[newRow];
		if (cEnd - cBegin == nBegin - nEnd) {
			for (int d = 0; d < cEnd - cBegin; d++) {
				temp = col_ind[row_ptr_begin[currentRow] + d];
				col_ind[row_ptr_begin[currentRow] + d] = col_ind[row_ptr_begin[newRow] + d];
				col_ind[row_ptr_begin[newRow] + d] = temp;
			}
			for (int d = 0; d < cEnd - cBegin; d++) {
				temp = values[row_ptr_begin[currentRow] + d];
				values[row_ptr_begin[currentRow] + d] = values[row_ptr_begin[newRow] + d];
				values[row_ptr_begin[newRow] + d] = temp;
			}
		}
		else {
			printf("difference in length");
			return false;
		}

	}
	return true;
}

static void offsetMatrix(int lastIndex, int lastRow, int row, int index, int upper,  int column, double value, bool search) {
	bool skip = false;
	if (search) {
		for (int s = index-1; s < upper && !skip; s++) {
			if (col_ind[s] < column && col_ind[s+1] > column) {
				column = col_ind[s]+1;
				index = s+1;
				skip = true;
				break;
			}
		}
	}
	if (index > row_ptr_end[row]) {
		index = row_ptr_end[row];
	}
	for (int l = lastIndex; l>= index; l--) {
		values[l + 1] = values[l];
	}
	for (int l = lastIndex; l>= index; l--) {
		col_ind[l + 1 ] = col_ind[l];
	}
	for (int r = row; r <= lastRow; r++ ) {
		row_ptr_end[r] += 1;
	}
	for (int r = row+1; r <= lastRow ; r++) {
		row_ptr_begin[r] += 1;
	}
	col_ind[index] = column;
	values[index] = value;
	
}
//Most important function
static bool pivot(int count, int countC) {
	int lastIndex = 0;
	//get length of values-array
	for (int i = 0; i < max_n_elements; i++) {
		if (values[i] != 0) {
			lastIndex = i;
		}
	}
	//get length of row_ptr array
	int lastRow = 0;
	for (int i = 0; i < max_n_rows; i++) {
		if (row_ptr_begin[i] != 0) {
			lastRow = i;
		}
	}
	//declare variables 
	double pivot, mult;
	//for each column, get pivot variable and perform gaussian elimination with it
	for (int i = 0; i < n_colsA; i++) {
		//get pivot value for current column
		pivot = values[row_ptr_beginU[i]+i];
		//go down from pivot value for each row
		for (int j = i + 1; j < n_rowsA; j++) {

			int col = col_ind[row_ptr_begin[j] + i];
			//If column of pivot variable not same as i, the value is 0 --> skip
			if (col != i) {
				continue;
			}
			//determine value to write into L-matrix and write it there 
			mult = values[row_ptr_begin[j] + i] / pivot;
			values[row_ptr_begin[j] + i] = mult;
			col_indLU[row_ptr_begin[j] + i] = i;
			//for each row j, perform A[j,k] = A[j,k]-mult*A[i,k]
			for (int k = j; k < n_colsA; k++) {
				//If k is outside of the matrix boundaries --> row finished for both i and j
				if (k > (col_ind[row_ptr_end[i]]- col_ind[row_ptr_begin[i]]) && k > (col_ind[row_ptr_end[j]]- col_ind[row_ptr_begin[j]])) {
					continue;
				}
				//k greater than i --> values[j] should stay the same
				if (k > col_ind[row_ptr_end[i]]- col_ind[row_ptr_begin[i]]) {
					continue;
				}
				//If k outside of j row boundaries
				if(row_ptr_begin[j]+k > row_ptr_end[j]){
					//Insert A[i,k] value to the beginning of A[j,k] and update following and dependent values 
					if (col_ind[row_ptr_begin[j]] > col_ind[row_ptr_begin[i]+k]) {
						offsetMatrix(lastIndex, lastRow, j, row_ptr_begin[j], 0, col_ind[row_ptr_begin[i] + k], values[row_ptr_begin[i] + k] * -1, false);
					}
					else {
						//insert A[i,k] value into jth row
						offsetMatrix(lastIndex, lastRow, j, row_ptr_begin[j], row_ptr_end[j], col_ind[row_ptr_begin[i] + k], values[row_ptr_begin[i] + k] * -1, true);
					}
					continue;
				}
				//If k indices of i and j different --> i or j starts earlier
				if (col_ind[row_ptr_begin[j] + k] != col_ind[row_ptr_begin[i] + k]) {
					int cdiff = col_ind[row_ptr_begin[i]] - col_ind[row_ptr_begin[j]];
					int cdiff2 = col_ind[row_ptr_end[i]] - col_ind[row_ptr_end[j]];
					//if i before j
					if (cdiff < 0) {
						offsetMatrix(lastIndex, lastRow, j, row_ptr_begin[j] + k, 0, col_ind[row_ptr_begin[i] + k], values[row_ptr_begin[i] + k ] *-1, false);
					}//else
					else {
						if (col_ind[row_ptr_begin[i] + k] < k) {
							continue;
						}
						else {
							values[row_ptr_begin[j] + k] = values[row_ptr_begin[j] + k] - mult * values[row_ptr_begin[j] + k];
						}
						
					}
				}
				else {
					//k values are the same for i and j
					values[row_ptr_begin[j] + k] = values[row_ptr_begin[j] + k] - mult * values[row_ptr_begin[i] + k];
					col_indU[row_ptr_beginU[j] + k] = k;
				}
			}

		}
	}
	return true;
	
}
//Split values matrix into two matrices: L for lower triangle and U for upper triangle
static void LUfactor() {
	//L matrix
	//init first row of L matrix
	row_ptr_beginL[0] = 0;
	row_ptr_endL[0] = row_ptr_beginL[0];
	//diagonal 1-matrix
	valuesL[0] = 1;
	col_indL[0] = 0;
	//continue for rest of L matrix, row ends one column "earlier" for each row
	for (int i = 1; i < n_rowsA; i++) {
		row_ptr_beginL[i] = row_ptr_endL[i-1] + 1;
		valuesL[row_ptr_beginL[i]+i] = 1;
		col_indL[row_ptr_beginL[i] + i] = i;
		row_ptr_endL[i] = row_ptr_beginL[i] + i;
		for (int j = 0; j < i &&  row_ptr_begin[i] + j <= row_ptr_end[i]; j++) {
			if (values[row_ptr_begin[i] + j] == 0)
				break;
			valuesL[row_ptr_beginL[i] + j] = values[row_ptr_begin[i] + j];
			col_indL[row_ptr_beginL[i] + j] = j;
		}
	}
	// U matrix
	//init first row of U matrix
	row_ptr_beginU[0] = 0;
	row_ptr_endU[0] = row_ptr_end[0];
	//continue for rest of U matrix, matrix starts one column "later" for each row
	for (int i = row_ptr_beginU[0]; i <= row_ptr_endU[0]; i++) {
		valuesU[i] = values[i];
		col_indU[i] = col_ind[i];
	}
	for (int i = 1; i < n_rowsA; i++) {
		row_ptr_beginU[i] = row_ptr_endU[i-1] +1;
		int offset = 0;
		int skip = 0;
		for (int k = row_ptr_begin[i]; k <= row_ptr_end[i]; k++) {
			if (col_ind[k] < i) {
				skip++;
			}
			else {
				break;
			}
		}
		for (int j = row_ptr_begin[i]+skip; j <= row_ptr_end[i]; j++) {
			
			valuesU[row_ptr_beginU[i]+offset] = values[j];
			col_indU[row_ptr_beginU[i]+offset] = offset;
			offset++;
		}
		row_ptr_endU[i] = row_ptr_beginU[i] + offset-1;
	}

}
static bool factorize(int count, int countC) {
	int arr;
	pBsize = sizeof(row_ptr_begin) / sizeof(*row_ptr_begin);
	arr = 10;
	//Switch rows so that no pivot value is zero
	int currentRow = -1, newRow = -1;
	currentRow, newRow = PivotIsZero();
	if (!switchRows(currentRow, newRow)) {
		printf("error");
		return false;
	}
	
	//compute A matrix from which L&U matrices are derived
	pivot(count, countC);
	//Initialize L & U matrices
	LUfactor();
	

	return true;

}
//Depending on index multiply different vectors and matrices, incomplete
static void multiplication(int index) {
	//L*(U*x) = b
	double prod;
	if (index == 0) {
		//for each row
		for (int i = 0; i < n_rowsA; i = i + 1)
		{
			//for each column
			for (int k = row_ptr_beginL[i]; k <= row_ptr_endL[i]; k = k + 1)
			{
				//for each column with a value in ValueL
				for (int c = 0; c < n_colsA; c++) {
					if (c == col_ind[k]) {
						prod = resultsx[c] * valuesL[row_ptr_begin[i] + k];
						resultsb[i] = resultsb[i] + prod;
					}
				}

			}
		}
	}
	//(U*X)
	if (index == 1) {
		//for each row
		for (int i = 0; i < n_rowsA; i = i + 1)
		{
			//for each column
			for (int k = row_ptr_beginU[i]; k <= row_ptr_endU[i]; k = k + 1)
			{
				//for each column with a value in ValueU
				for (int c = 0; c < n_colsA; c++) {
					if (c == col_indU[k]) {
						prod = xvector[c] * valuesU[row_ptr_begin[i] + k];
						resultsx[i] = resultsx[i] + prod;
					}
				}
			}
		}
	}
	//A*x
	if (index == 2) {
		//for each row
		for (int i = 0; i < n_rowsA; i = i + 1)
		{
			//for each column
			for (int k = row_ptr_begin[i]; k <= row_ptr_end[i]; k = k + 1)
			{
				//for each column with a value in ValueU
				for (int c = 0; c < n_colsA; c++) {
					if (c == col_ind[k]) {
						prod = xvector[c] * values[row_ptr_begin[i] + k];
						results[i] = results[i] + prod;
					}
				}
			}
		}
	}
}
//No gaussian substitution implemented
static void substitution(int index) {
	//index = 0 --> Ly = resultB
	int col_in = 9;
}

//Initialize x-vector
//index = 0 --> all 1
//index = 1 --> all .1
//index = 2 --> alternating 1, -1
//index = 3 --> alternating 5, -5
//index = 4 --> alternating 100, -100
static void init(int index) {
	if (index == 0) {
		for (int m = 0; m < sizeof(xvector) / sizeof(xvector[0]); m++) {
			xvector[m] = 1;
			xvector2[m] = 1;
		}
	}
	if (index == 1) {
		for (int m = 0; m < sizeof(xvector) / sizeof(xvector[0]); m++) {
			xvector[m] = .1;
			xvector2[m] = .1;
		}
	}
	if (index == 2) {
		for (int m = 0; m < sizeof(xvector) / sizeof(xvector[0]); m += 2) {
			xvector[m] = 1;
			xvector[m + 1] = -1;
			xvector2[m] = 1;
			xvector2[m + 1] = -1;
		}
	}
	if (index == 3) {
		for (int m = 0; m < sizeof(xvector) / sizeof(xvector[0]); m += 2) {
			xvector[m] = 5;
			xvector[m + 1] = -5;
			xvector2[m] = 5;
			xvector2[m + 1] = -5;
		}
	}
	if (index == 4) {
		for (int m = 0; m < sizeof(xvector) / sizeof(xvector[0]); m += 2) {
			xvector[m] = 100;
			xvector[m + 1] = -100;
			xvector2[m] = 100;
			xvector2[m + 1] = -100;
		}
	}

}
/*

int
main(int argc, char **argv)
{
	int nnz, n_rows, n_cols;
	bool ok(false);
	//change argv[1] to filename if you wish to run the script from an IDE with a fixed matrix
	const char * filename = "mcfe.mtx";
	ok = load_matrix_market(filename, max_n_elements, max_n_rows,
		nnz, n_rows, n_cols,
		values, col_ind, row_ptr_begin, row_ptr_end);
	n_rowsA = n_rows;
	n_colsA = n_cols;
	for (int i = 0; i < max_n_elements; i++) {
		valuesLU[i] = values[i];
		col_indLU[i] = col_ind[i];
	}
	for (int i = 0; i < max_n_rows; i++) {
		row_ptr_beginLU[i] = row_ptr_begin[i];
		row_ptr_endLU[i] = row_ptr_end[i];
	}

	
	if (!ok)
	{
		fprintf(stderr, "failed to load matrix.\n");
		return -1;
	}

	for (int i = 1; i < n_cols; i++) {
		int maxCol = 0;
		int maxColIndex = -1;
		int  maxRowIndex = -1;
		int sum = 0;
		for (int r = 1; r < pBsize; r++) {
			if (maxCol < values[row_ptr_begin[r] + i - 1]) {
				maxCol = values[row_ptr_begin[r] + i - 1];
				maxColIndex = row_ptr_begin[r] + i - 1;
				maxRowIndex = r + i - 1;
			}
			//sparse matrix L put values of best row to the top

			//sparse matrix U put values of best row to the bottom
		}
	}


	int COUNTV = 0;
	int Count2 = 0;
	int vSize = sizeof(values) / sizeof(values[0]);
	for (int i = 0; i < vSize; i++) {
		if (values[i] != 0) //check if element non-zero
		{
			COUNTV++;
		}
		if (values[i] == 0) {
			Count2++;
		}
	}
	int COUNTC = 0;
	int cSize = sizeof(col_ind) / sizeof(col_ind[0]);
	for (int i = 0; i < cSize; i++) {
		if (col_ind[i] != 0) // check if element non - zero
			COUNTC++;
	}
	//initialize different x vectors
	init(0);


	/* For debugging, can be removed when implementation is finished.*/
	//dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);
	//bool(stop) = false;
/*
	struct timespec start_time;
	clock_gettime(CLOCK_REALTIME, &start_time);

	/* Perform LU factorization here
	//check if any errors occurred
	if (!factorize(COUNTV, COUNTC)) {
		exit(-2);
	}
	//call different multiplication executions, according to index 2 => A*x = b
	multiplication(2);
	//No forward/backward substitution yet


	struct timespec end_time;
	clock_gettime(CLOCK_REALTIME, &end_time);


	struct timespec elapsed_time;
	timespec_subtract(&elapsed_time, &end_time, &start_time);

	double elapsed = (double)elapsed_time.tv_sec +
		(double)elapsed_time.tv_nsec / 1000000000.0;
	fprintf(stderr, "elapsed time: %f s\n", elapsed);


	return 0;
}*/