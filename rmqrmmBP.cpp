//============================================================================
// Name        : rmqrmmBP.cpp
// Author      : Hector Ferrada
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "includes/RMQRMM64.h"

using namespace std;
using namespace rmqrmm;

bool TRACE = false;		// true: print all details for console
bool TEST = true;		// true: apply exhaustive test
uint N_REP = 100000;		// number of REPETITIONS
bool RANDOM = true;		// true, random vales (i,j) for test rmq->queryRMQ(i,j)
bool DOWN = true;		// true, decrease pseudosorted arrays
bool LOAD_CREATE = true;		// 0/1 LOAD/CRATE

int testRMQ(int argc, char *argv[]);
int testRMQLoad(int argc, char *argv[]);
double runExperimentTimeRMQ(RMQRMM64 *rmq, ulong n);

int main(int argc, char *argv[]) {

	if(argc == 2){
		testRMQLoad(argc, argv);
	}else
		testRMQ(argc, argv);
}

int testRMQ(int argc, char *argv[]) {
	RMQRMM64 *rmq = NULL;
	uint pseudoSorted;
	char resultsFile[400], saveLoadFile[400];
	long int *A;
	ulong n, i;
	double ct_rmq, qt_rmq;

	ct_rmq = 0.0;

	if(argc != 8){
		cout << "Missing parameters..." << endl;
		cout << "Fischer DFUDS-RMQ's usage: " << endl;
		cout << " " << argv[0] << " <array length(n)> <load/create 0/1> <saveLoadFile> <repetitions> <pseudoSorted> <RandomWeight> <resultsFile>" << endl;
		return 1;
	}

	n = atol(argv[1]);
	LOAD_CREATE = atoi(argv[2]); // 0/1
	strcpy(saveLoadFile, argv[3]);
	N_REP = atol(argv[4]);
	pseudoSorted = atoi(argv[5]);
	ulong weight = atol(argv[6]);
	strcpy(resultsFile, argv[7]);

	cout << "Parameters.." << endl;
	cout << " LOAD_CREATE = " << LOAD_CREATE << endl;
	cout << " saveLoadFile = " << saveLoadFile << endl;
	cout << " N_REP = " << N_REP << endl;
	cout << " pseudoSorted = " << pseudoSorted << endl;
	cout << " weight = " << weight << endl;
	cout << " resultsFile = " << resultsFile << endl;
	cout << " pseudoSorted DOWN = " << DOWN << endl;

	if(LOAD_CREATE){ // create...
		A = new long int[n];
		cout << "Creating the sequence A[1.." << n << "] with numbers in the range [0.." << weight-1 << "]..." << endl;

		if(pseudoSorted){
			if (DOWN){
				// decrease
				for(size_t i=0; i < n; ++i)
					A[i] = n-i-weight + rand()%(2*weight);
			}else{
				// increase
				for(size_t i=0; i < n; ++i)
					A[i] = (i-weight + rand()%(2*weight));
			}
		}else{
			for(size_t i=1; i < n; ++i)
				A[i] = (rand()%weight);
		}

		{
			// Saving A[] for test in 'fileName'
			char fileName[400];
			strcpy(fileName, saveLoadFile);
			strcat(fileName, ".Seq");
			cout << "Saving A[] for test in " << fileName << endl;
			ofstream os (fileName, ios::binary);
			os.write((const char*)A, n*sizeof(long int));
			os.close();
		}

		if (TRACE){
			cout << "A[0.." << n-1 << "]" << endl;
			for (i=0; i<n; i++){
				cout << A[i] << " ";
				if ((i+1)%20==0)
					cout << endl;
			}
		}

		cout << "Creating rmq for A[]..." << endl;
		ct_rmq = getTime_ms();
		rmq = new RMQRMM64(A, n);
		ct_rmq = getTime_ms() - ct_rmq;
		cout << "rmq creted in " << ct_rmq << " milliseconds" << endl;
		rmq->saveDS(saveLoadFile);
		rmq->~RMQRMM64();
		rmq->loadDS(saveLoadFile);
	}else	// load...
		rmq = new RMQRMM64(saveLoadFile);

	cout << "rmq size = " << rmq->getSize() << " Bytes = " << (float)rmq->getSize()*8.0/(float)n << " n" << endl;

	if (TEST){
		cout << "Test RMQ..." << endl;

		if(LOAD_CREATE == 0)		{
			// Loading A[] for test in 'fileName'
			A = new long int[n];
			char fileName[400];
			strcpy(fileName, saveLoadFile);
			strcat(fileName, ".Seq");
			cout << "Loading A[] for test in " << fileName << endl;
			ifstream is(fileName, ios::binary);
			is.read((char*)A, n*sizeof(long int));
			is.close();
		}

		ulong min, rmq_min, beg, end;

		/*beg=92;
		end=924;
		for (min=beg, i=beg+1; i<=end; i++){
			if (A[i] < A[min])
				min = i;
		}
		rmq_min = rmq->queryRMQ(beg,end);
		if (rmq_min < beg || rmq_min > end){
			cout << "ERROR... rmq_min = " << rmq_min << " out of range [" << beg << " , " << end << " ]" << endl;
			exit(1);
		}else{
			if (A[rmq_min] != A[min]){
				cout << "ERROR... rmq(" << beg << " , " << end << " ) = " << rmq_min << " != " << min << endl;
				exit(1);
			}
		}*/

		if (RANDOM){
			cout << "Random test O(N_REP)..." << endl;
			for (uint t=0; t<N_REP; t++){
				beg = (rand() % (n/2));
				end = n/2 + (rand() % (n/2)-1);
				if (end > beg+1){
					for (min=beg, i=beg+1; i<=end; i++){
						if (A[i] < A[min])
							min = i;
					}
					if (beg==end)
						rmq_min = beg;
					else
						rmq_min = rmq->queryRMQ(beg,end);
					if (rmq_min < beg || rmq_min > end){
						cout << "ERROR... rmq_min = " << rmq_min << " out of range [" << beg << " , " << end << " ]" << endl;
						exit(1);
					}else{
						if (A[rmq_min] != A[min]){
							cout << "ERROR... (" << beg << " , " << end << " ) = " << rmq_min << " != " << min << endl;
							exit(1);
						}
					}
				}
			}
		}else{
			cout << "Exhaustive test O(n^2)..." << endl;
			for (beg=0; beg<n; beg++){
				min = beg;
				for (end=beg; end<n; end++){
					if (A[end] < A[min])
						min = end;
					rmq_min = rmq->queryRMQ(beg,end);

					if (rmq_min < beg || rmq_min > end){
						cout << "ERROR... rmq_min = " << rmq_min << " out of range [" << beg << " , " << end << " ]" << endl;
						exit(1);
					}else{
						if (A[rmq_min] != A[min]){
							cout << "ERROR... (" << beg << " , " << end << " ) = " << rmq_min << " != " << min << endl;
							exit(1);
						}
					}
				}
			}
		}

		cout << "rmq test completed!" << endl;
	}

	if(LOAD_CREATE) // crate..
		delete[] A;

	qt_rmq = runExperimentTimeRMQ(rmq, n);

	FILE *fp = fopen(resultsFile, "a+" );
	// [n] [weight] [o(n) extra size] [construction time ms] [query time ms]
	fprintf(fp, "%lu %lu %f %G %G\n", n, weight, (float)rmq->getSize()*8.0/(float)n, ct_rmq, qt_rmq);
	fclose(fp);

	rmq->~RMQRMM64();
	cout << "Program Finished OK" << endl;	return 0;
}


int testRMQLoad(int argc, char *argv[]) {
	cout << "Loading..." << endl;
	char* file = new char[100];
	strcpy(file, argv[1]);
	RMQRMM64 *rmq = new RMQRMM64(file);

	cout << "rmq size = " << rmq->getSize() << " Bytes = " << (float)rmq->getSize()*8.0/((float)rmq->nP/2.0) << " n" << endl;
	runExperimentTimeRMQ(rmq, rmq->nP/2);
	rmq->~RMQRMM64();
	cout << "Program Finished OK" << endl;	return 0;
}

double runExperimentTimeRMQ(RMQRMM64 *rmq, ulong n){
	double t, avgTime;
	ulong i, j;
	uint k;

	cout << "_______________________________________________________________________________________________________" << endl;
	cout << "Start Experiment (random values for rmq(i,j))..." << endl;
	avgTime = 0.0;

	for (k=0; k<N_REP; k++){
		i = (rand() % (n-1));
		j = i + (rand() % (n-(i+1)));
		t = getTime_ms();
		rmq->queryRMQ(i,j);
		t = getTime_ms() - t;
		avgTime += t;
	}

	avgTime /= (double)N_REP;
	cout << "Average CPU time for item : " << avgTime*1000.0 << " Microseconds" << endl;
	cout << "_______________________________________________________________________________________________________" << endl;

	return avgTime;
}
