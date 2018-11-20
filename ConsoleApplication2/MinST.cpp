#include "stdafx.h"
#include "matrix.h"
 #include <tuple>
#include <set>
#include <algorithm>
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



struct Edge {
	int src;
	int target;
	float weigth;
};

struct MinST {
	std::vector<Edge*> edges;
};

struct Neighblist {
	std::set<int> parentNodes;
	std::vector<int> nodelist;
	std::vector<float> weigthlist;
};
struct NodeList
{
	std::vector<Neighblist*> nodes;
};
static NodeList* nodelist;
static NodeList* nodelistOrig;
static vector<std::tuple<int, int>> pairs;
static vector<std::tuple<int, int>> MST;
static MinST Tree;
static int pairsize;
static int nodelistsize;

NodeList* toNodeList(int n_rows, int n_cols) {
	NodeList* nodes = new NodeList;
	for (int i = 0; i < sizeof(row_ptr_begin); i++) {
		Neighblist* nList = new Neighblist;
		for (int j = 0; j <= row_ptr_end[i] - row_ptr_begin[i] && row_ptr_end[i] != 0; j++) {
			if (col_ind[row_ptr_begin[i] + j] == i)
				continue;
			nList->parentNodes.insert(i);
			nList->nodelist.push_back(col_ind[row_ptr_begin[i] + j]);
			nList->weigthlist.push_back(values[row_ptr_begin[i] + j]);
		}
		int len = (nList->nodelist.size());
		if (nList->nodelist.size() > 0) {
			nodes->nodes.push_back(nList);
		}
	}
	return nodes;
}

int findMinEdge(int j, Neighblist* nodes, vector<int> contractedlist) {
	float temp = 100000;
	int index = -1;
	for (int i = 0; i < nodes->weigthlist.size(); i++) {
		auto test = std::find(contractedlist.begin(), contractedlist.end(), int(nodes->nodelist.at(i)));
		bool found = test == contractedlist.end();
		if (nodes->weigthlist.at(i) < temp && found) {
			temp = nodes->weigthlist.at(i);
			index = nodes->nodelist.at(i);
		}
	}
	if (index < 0)
		index = j;
	return index;
}

float findWeight(int src, int target) {
	Neighblist* node1 = nodelist->nodes.at(src);
	int index = 0;
	for (int i = 0; i < node1->nodelist.size(); i++) {
		if (node1->nodelist.at(i) == target) {
			index = i;
			break;
		}
	}
	return node1->weigthlist.at(index);
}

void getPairs() {
	pairs.clear();
	std::vector<int> contractedList;
	for (int i = 0; i < nodelist->nodes.size(); i++) {
		if (i == 92)
			char stop = 'sfg';
		if (std::count(contractedList.begin(), contractedList.end(), i) != 0)
			continue;
		int index = findMinEdge(i, nodelist->nodes.at(i), contractedList);
		if (std::count(contractedList.begin(), contractedList.end(), index) != 0)
			continue;
		contractedList.push_back(i);
		contractedList.push_back(index);
		pairs.push_back(std::make_tuple(i, index));
		Edge* e = new Edge;
		e->src = i;
		e->target = index;
		e->weigth = findWeight(i, index);
		Tree.edges.push_back(e);
		string s = "m";
	}
	pairsize = pairs.size();
	//findMinEdge for each neighblist
	//create pairs
	//merge pairs
	//save as nodelist
}

int lookup(int node) {
	for (int i = 0; i < pairsize; i++) {
		std::tuple<int, int> tp = pairs[i];
		if (std::get<0>(tp) == node) {
			return node;
		}
		else if (std::get<1>(tp) == node) {
			return std::get<0>(tp);
		}
	}
}
Neighblist* findNode(int index) {
	for (int i = 0; i < nodelist->nodes.size(); i++) {
		std::set<int> temp = nodelist->nodes[i]->parentNodes;
		if (temp.find(index) != temp.end()) {
			return nodelist->nodes[i];
		}
	}
}
void Union() {
	NodeList* nodes2 = new NodeList;
	for (int i = 0; i < pairsize; i++) {
		std::tuple<int, int> pair = pairs[i];
		Neighblist* node1 = findNode(std::get<0>(pair));
		Neighblist* node2 = findNode(std::get<1>(pair));
		if (node1 == node2) {
			nodes2->nodes.push_back(node1);
			continue;
		}
		Neighblist* nList = new Neighblist;
		
		nList->parentNodes.insert(node1->parentNodes.begin(), node1->parentNodes.end());
		nList->parentNodes.insert(node2->parentNodes.begin(), node2->parentNodes.end());
						
		int n1size = node1->nodelist.size();
		int n2size = node2->nodelist.size();
		for (int j = 0; j < n1size; j++) {
			int l1 = lookup(node1->nodelist[j]);
			if (node1->nodelist[j] == std::get<1>(pair)) {
				node1->nodelist.erase(node1->nodelist.begin() + j);
				node1->weigthlist.erase(node1->weigthlist.begin() + j);
				n1size--;
			}
			if (j< n1size && l1 != node1->nodelist.at(j) && std::find(node1->nodelist.begin(), node1->nodelist.end(), l1) != node1->nodelist.end()){
				int index = std::find(node1->nodelist.begin(), node1->nodelist.end(), l1)- node1->nodelist.begin();
				float weight = min(node1->weigthlist.at(j), node1->weigthlist.at(index));
				node1->weigthlist.at(index) = weight;
				node1->nodelist.erase(node1->nodelist.begin() + j);
				node1->weigthlist.erase(node1->weigthlist.begin() + j);
				n1size--;
				j--;
			}
		}
		for (int j = 0; j < n2size; j++) {
			int l2 = lookup(node2->nodelist.at(j));
			if (node2->nodelist.at(j) == std::get<0>(pair)) {
				node2->nodelist.erase(node2->nodelist.begin() + j);
				node2->weigthlist.erase(node2->weigthlist.begin() + j);
				n2size--;
			}
			if (j< n2size && l2 != node2->nodelist.at(j) && std::find(node2->nodelist.begin(), node2->nodelist.end(), l2) != node2->nodelist.end()) {
				int index = std::find(node2->nodelist.begin(), node2->nodelist.end(), l2) != node2->nodelist.end();
				float weight = min(node2->weigthlist.at(j), node2->weigthlist.at(index));
				node2->weigthlist.at(index) = weight;
				node2->nodelist.erase(node2->nodelist.begin() + j);
				node2->weigthlist.erase(node2->weigthlist.begin() + j);
				n2size--;
				j--;
			}
		}
		//merge lists
		int a = 0;
		int b = 0;
		int node1add = 0;
		int node2add = 0;
		while ((a < (n1size) && b < (n2size))) {
			int la = lookup(node1->nodelist.at(a));
			int lb = lookup(node2->nodelist.at(b));
			if (node1add < n1size && (la < lb || (b == (n2size-1) && la != lb))) {
				nList->nodelist.push_back(la);
				nList->weigthlist.push_back(node1->weigthlist.at(a));
				a < (n1size-1) ? a++ : a;
				node1add++;
			}

			else if (node2add < n2size && (la > lb || (a == (n1size - 1) && la != lb))) {
				nList->nodelist.push_back(lb);
				nList->weigthlist.push_back(node2->weigthlist.at(b));
				b < (n2size-1) ? b++ : b;
				node2add++;
			}
			else{
				if (node1->weigthlist.at(a) < node2->weigthlist.at(b)) {
					nList->nodelist.push_back(la);
					nList->weigthlist.push_back(node1->weigthlist.at(a));

				}
				else {
					nList->nodelist.push_back(lb);
					nList->weigthlist.push_back(node2->weigthlist.at(b));
				}
				a < (n1size - 1) ? a++ : a;
				b < (n2size - 1) ? b++ : b;
				node1add++;
				node2add++;
			}
			if (node1add >= n1size  && node2add >= n2size) {
				break;
			}
		}
		//if (nodes2->nodes.size() == 1100) {
		
		nodes2->nodes.push_back(nList);
	}
	nodelist = nodes2;
}


int main(int argc, char **argv)
{
	int nnz, n_rows, n_cols;
	bool ok(false);
	//change argv[1] to filename if you wish to run the script from an IDE with a fixed matrix
	const char * filename = "nopoly.mtx";
	ok = load_matrix_market(filename, max_n_elements, max_n_rows,
		nnz, n_rows, n_cols,
		values, col_ind, row_ptr_begin, row_ptr_end);

	//Convert matrix to struct
	nodelist = toNodeList(n_rows, n_cols);
	nodelistsize = nodelist->nodes.size();
	int nodelistsizeOld = 0;
	//merge pairs
	//Find number of components
	while(nodelistsize > 1000 || nodelistsize != nodelistsizeOld){
	nodelistsize = nodelist->nodes.size();
	getPairs();
	Union();
	int printS = 'S';
	nodelistsizeOld = nodelist->nodes.size();
	}
	/* For debugging, can be removed when implementation is finished.*/
	//dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);
	//bool(stop) = false;

	struct timespec start_time;
	clock_gettime(CLOCK_REALTIME, &start_time);



	struct timespec end_time;
	clock_gettime(CLOCK_REALTIME, &end_time);


	struct timespec elapsed_time;
	timespec_subtract(&elapsed_time, &end_time, &start_time);

	double elapsed = (double)elapsed_time.tv_sec +
		(double)elapsed_time.tv_nsec / 1000000000.0;
	fprintf(stderr, "elapsed time: %f s\n", elapsed);


	return 0;
}