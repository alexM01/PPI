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
static int row_ptr_begin[max_n_rows];
static int row_ptr_end[max_n_rows];
static int pBsize;
static double valuesL[131072];
static double valuesU[131072];
static double col_indL[131072];
static double col_indU[131072];
static int row_ptr_beginL[16384];
static int row_ptr_beginU[16384];
static int  row_ptr_endL[16384];
static int row_ptr_endU[16384];
static int skip[16384];
static int skip2[16384];
static double results[131072];
static double resultsx[131072];
static double resultsb[131072];
static double xvector[131072];
static double xvector2[131072];

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

static void factorL() {
	bool newline = true;
	double pivot, mult;
	for (int i = 0; i < n_colsA; i++) {
		pivot = values[row_ptr_begin[i] + i];
		for (int j = i + 1; j < n_rowsA; j++) {
			int stop3, stop4;
			if (row_ptr_beginL[j] + i == 155) {
				stop3 = row_ptr_beginL[j] + i;
				stop4 = row_ptr_endL[j];
				int stop5 = 9;
			}
			if ((row_ptr_beginL[j] + i) > 0 && (row_ptr_beginL[j] + i) >= row_ptr_endL[j]) {
				int diff = (row_ptr_endL[j] - row_ptr_beginL[j]);
				continue;
			}
			if ((i > col_ind[row_ptr_begin[j] + i]) || (i > (row_ptr_end[j] - row_ptr_begin[j]))) {
				int diff = (row_ptr_end[j] - row_ptr_begin[j]);
				continue;
			}
			if (values[row_ptr_begin[j] + i] == 0) {
				skip2[j++];
				continue;
			}
			mult = values[row_ptr_begin[j] + i] / pivot;
			if (i == 1 && j == 17) {
				int stop = values[row_ptr_begin[j] + i];
			}
			if (newline && mult != 0 && i == 0) {


				row_ptr_beginL[j] = row_ptr_endL[j - 1] + 1;
				//save end pointer number, subtract 1 for each 0 per row
				row_ptr_endL[j] = row_ptr_beginL[j] + j - skip[j];
				valuesL[row_ptr_endL[j]] = 1;

				col_indL[row_ptr_endL[j]] = j - skip2[j];
			}

			if (j > skip[j]) {
				valuesL[row_ptr_beginL[j] + i] = mult;
				//valuesU[row_ptr_begin[j] + i] = 0;
				col_indL[row_ptr_beginL[j] + i] = i + skip[j];
			}

		}
	}
}

static bool pivot(int count, int countC) {

	valuesL[0] = 1;
	col_indL[0] = 0;
	row_ptr_beginL[0] = 0;
	row_ptr_endL[0] = 0;
	//copy first row
	for (int k = 0; k <= row_ptr_end[0]; k++) {
		valuesU[k] = values[k];
		col_indU[k] = col_ind[k];
		if (k > 0) {
			valuesL[k] = 0;
		}
	}
	row_ptr_beginU[0] = row_ptr_begin[0];
	row_ptr_endU[0] = row_ptr_end[0];
	double pivot, mult;
	bool newline = true;
	bool red = false;
	int skipK = 0;
	for (int i = 0; i < n_colsA; i++) {
		pivot = values[row_ptr_begin[i] + i];
		int skipu = 0;

		for (int j = i + 1; j < n_rowsA; j++) {
			if (i > (row_ptr_end[j] - row_ptr_begin[j])) {
				int diff = (row_ptr_end[j] - row_ptr_begin[j]);
				skipK += row_ptr_end[j] - row_ptr_begin[j];
				continue;
			}

			mult = values[row_ptr_begin[j] + i] / pivot;
			if (i == 0 && col_ind[row_ptr_begin[j]] > 0) {
				skip[j] += j > col_ind[row_ptr_begin[j]] ? col_ind[row_ptr_begin[j]] : j;
				
			}

			if (col_ind[row_ptr_begin[j] + i] < col_ind[row_ptr_begin[i] + i]) {
				int sdiff = col_ind[row_ptr_begin[j] + i] - col_ind[row_ptr_begin[i] + i];
				skip[j] += sdiff;
				continue;
			}
			else {

			}
			int offset = 0;
			row_ptr_beginU[j] = row_ptr_endU[j - 1] + 1;
			if (j ==16) {
				int s = 1;
			}

			for (int k = j; k < n_colsA; k++) {

				int diffJ = row_ptr_end[j] - row_ptr_begin[j];
				int diffI = row_ptr_end[i] - row_ptr_begin[i];
				int maxDiff = diffJ - diffI > 0 ? diffJ : diffI;
				
				if (k > maxDiff && row_ptr_begin[j]+k > row_ptr_end[j] && col_ind[row_ptr_begin[j]] < k && col_ind[row_ptr_end[j]] > k) {
					int diffK = k - maxDiff;
					continue;
				}
				//row_ptr_beginU[j] = row_ptr_endU[j - 1] + 1;

				if (col_ind[row_ptr_begin[j] + k ] != col_ind[row_ptr_begin[i] + k]) {
					int diff1 = col_ind[row_ptr_begin[j + 1] + k];
					int diff2 = col_ind[row_ptr_begin[j] + k];
					int cdiff = col_ind[row_ptr_begin[j + 1] + k] - col_ind[row_ptr_begin[j] + k];
					if (col_ind[row_ptr_begin[j]] < k) {
						row_ptr_endU[j] = row_ptr_beginU[j] + (row_ptr_end[j] - row_ptr_begin[j] - k);
					}
					else {
						row_ptr_endU[j] = row_ptr_beginU[j] + (row_ptr_end[j] - row_ptr_begin[j]);
					}
						
					if (cdiff < 0) {
						valuesU[row_ptr_begin[j] + k] = values[row_ptr_begin[i] + k] * -1;
						col_indU[row_ptr_begin[j] + k] = col_ind[row_ptr_begin[i] + k];
						row_ptr_endU[j] = row_ptr_endU[j] + 1;
					}
					else {
						if (row_ptr_begin[j] + k < row_ptr_end[j]) {
							valuesU[row_ptr_begin[j] + k] = values[row_ptr_begin[j] + k];
							col_indU[row_ptr_begin[j] + k] = col_ind[row_ptr_begin[j] + k];
							row_ptr_endU[j] = row_ptr_endU[j] + 1;
						}
						else {
							for (int cd = 0; cd < abs(cdiff); cd++) {
								valuesU[row_ptr_begin[j] + cd] = values[row_ptr_begin[j] + k + cd];
								col_indU[row_ptr_begin[j] + cd] = col_ind[row_ptr_begin[j] + k + cd];
								row_ptr_endU[j] = row_ptr_endU[j] + 1;
							}
						}
					}
				}
				else {

					valuesU[row_ptr_beginU[j] + offset - skipK] = values[row_ptr_begin[j] + k] - mult * values[row_ptr_begin[i] + k];
					col_indU[row_ptr_beginU[j] + offset - skipK] = k;
					if (j == k) {
						if (valuesU[row_ptr_beginU[j] + offset - skipK] != 0) {
							//row_ptr_beginU[j] = row_ptr_endU[j - 1] + 1;
							row_ptr_endU[j] = row_ptr_beginU[j] + (row_ptr_end[j] - row_ptr_begin[j] - k - skipK);
						}
					}
					if (valuesU[row_ptr_beginU[j] + offset - skipK] == 0) {
						row_ptr_endU[j]--;
					}
				}

				offset++;
			}
			offset = 0;
			newline = true;
		}
	}
	factorL();
	return true;
}

static bool factorize(int count, int countC) {
	int arr;
	pBsize = sizeof(row_ptr_begin) / sizeof(*row_ptr_begin);
	arr = 10;/*
	int currentRow = -1, newRow = -1;
	currentRow, newRow = PivotIsZero();
	if (!switchRows(currentRow, newRow)) {
	printf("error");
	return false;
	}*/
	//Factorize
	//Initialize L & U matrices
	pivot(count, countC);
	/*
	int COUNT = 0;
	int vSize = sizeof(valuesL) / sizeof(valuesL[0]);
	for (int i = 0; i < vSize; i++) {
	if (valuesL[i] != 0) //which means this element is true, not empty element.
	COUNT++;
	}
	int COUNT1 = 0;
	int vSize1 = sizeof(valuesU) / sizeof(valuesU[0]);
	for (int i = 0; i < vSize1; i++) {
	if (valuesU[i] != 0) //which means this element is true, not empty element.
	COUNT1++;
	}*/

	return true;

}
static void multiplication(int index) {
	double prod;
	if (index == 0) {
		for (int i = 0; i < n_rowsA; i = i + 1)
		{
			for (int k = row_ptr_beginL[i]; k <= row_ptr_endL[i]; k = k + 1)
			{
				for (int c = 0; c < n_colsA; c++) {
					if (c == col_indL[k]) {
						prod = resultsx[c] * valuesL[row_ptr_begin[i] + k];
						resultsb[i] = resultsb[i] + prod;
					}
				}

			}
		}
	}

	if (index == 1) {
		for (int i = 0; i < n_rowsA; i = i + 1)
		{
			for (int k = row_ptr_beginU[i]; k <= row_ptr_endU[i]; k = k + 1)
			{
				for (int c = 0; c < n_colsA; c++) {
					if (c == col_indU[k]) {
						prod = xvector[c] * valuesU[row_ptr_begin[i] + k];
						resultsx[i] = resultsx[i] + prod;
					}
				}
			}
		}
	}
	if (index == 2) {
		for (int i = 0; i < n_rowsA; i = i + 1)
		{
			for (int k = row_ptr_begin[i]; k <= row_ptr_end[i]; k = k + 1)
			{
				for (int c = 0; c < n_colsA; c++) {
					if (c == col_ind[k]) {
						prod = xvector[c] * values[row_ptr_begin[i] + k];
						results[i] = results[i] + prod;
					}
				}
			}
		}
	}
}/*
static void MatrixMultiplication() {
for (int i = 0; i < n_rowsA; i = i + 1)
{
for (int k = row_ptr_beginL[i]; k <= row_ptr_endL[i]; k = k + 1)
{
for (int c = 0; c < n_colsA; c++) {
if (col_indU[k] == c) {
resultsM[i] = resultsM[i] + valuesU[k] * valuesL[row_ptr_begin[i]];
}

}

}
}
}*/

static void substitution(int index) {
	//index = 0 --> Ly = resultB
	int col_in = 9;
}


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


int
main(int argc, char **argv)
{
	printf("Tset");
	/*
	if (argc != 2)
	{
	fprintf(stderr, "usage: %s <filename>\n", argv[0]);
	return -1;
	}
	*/
	int nnz, n_rows, n_cols;
	bool ok(false);
	const char * filename = "mcfe.mtx";
	ok = load_matrix_market(filename, max_n_elements, max_n_rows,
		nnz, n_rows, n_cols,
		values, col_ind, row_ptr_begin, row_ptr_end);
	n_rowsA = n_rows;
	n_colsA = n_cols;
	if (!ok)
	{
		fprintf(stderr, "failed to load matrix.\n");
		return -1;
	}
	/*
	double** l = new double*[n_rows];
	for (int i = 1; i < n_rows; ++i)
	l[i] = new double[n_cols];


	double** u = new double*[n_rows];
	for (int i = 1; i < n_rows; ++i)
	u[i] = new double[n_cols];*/

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
		if (values[i] != 0) //which means this element is true, not empty element.
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
		if (col_ind[i] != 0) //which means this element is true, not empty element.
			COUNTC++;
	}
	init(0);

	//double nvalues[COUNT];



	//CSR Multiplication


	/*
	int COUNT = 0;
	int vSize2 = sizeof(valuesL) / sizeof(valuesL[0]);
	for (int i = 0; i < vSize2; i++) {
	if (valuesL[i] == 1) //which means this element is true, not empty element.
	COUNT++;
	}
	int COUNT1 = 0;
	int vSize1 = sizeof(valuesU) / sizeof(valuesU[0]);
	for (int i = 0; i < vSize1; i++) {
	if (valuesU[i] != 0) //which means this element is true, not empty element.
	COUNT1++;
	}int COUNTcL = 0;
	int cSizeL = sizeof(col_indL) / sizeof(col_indL[0]);
	for (int i = 0; i < cSizeL; i++) {
	if (col_indL[i] != 0) //which means this element is true, not empty element.
	COUNTcL++;
	}
	int COUNT1cU = 0;
	int cSizeU = sizeof(col_indU) / sizeof(col_indU[0]);
	for (int i = 0; i < cSizeU; i++) {
	if (col_indU[i] != 0) //which means this element is true, not empty element.
	COUNT1cU++;
	}
	int s = 1;*/

	/*
	Factorization
	int n, i, k, j, p, sum;
	for(k=1;k<=n;k++)
	{
	u[k][k]=1;
	for(i=k;i<=n;i++)
	{
	sum=0;
	for(p=1;p<=k-1;p++)
	sum+=l[i][p]*u[p][k];
	l[i][k]=a[i][k]-sum;
	}

	for(j=k+1;j<=n;j++)
	{
	sum=0;
	for(p=1;p<=k-1;p++)
	sum+=l[k][p]*u[p][j];
	u[k][j]=(a[k][j]-sum)/l[k][k];
	}
	}*/


	/* For debugging, can be removed when implementation is finished.*/
	//dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);
	//bool(stop) = false;

	struct timespec start_time;
	clock_gettime(CLOCK_REALTIME, &start_time);

	/* Perform LU factorization here*/

	if (!factorize(COUNTV, COUNTC)) {
		exit(-2);
	}
	//pb
	multiplication(1);
	//substitution

	int stop2 = 9;
	//MatrixMultiplication();
	int stop = 10;


	struct timespec end_time;
	clock_gettime(CLOCK_REALTIME, &end_time);


	struct timespec elapsed_time;
	timespec_subtract(&elapsed_time, &end_time, &start_time);

	double elapsed = (double)elapsed_time.tv_sec +
		(double)elapsed_time.tv_nsec / 1000000000.0;
	fprintf(stderr, "elapsed time: %f s\n", elapsed);


	return 0;
}