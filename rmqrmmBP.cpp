//============================================================================
// Name        : rmqrmmBP.cpp
// Author      : Hector Ferrada
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
// You can run this binary after to compile this with "make test", then you can call it with:
//	./rmqrmmBP 1000000 1 pathFileStructure.rmq 10000 0 0 10000 fileResults.txt pathResult/ 0
//

#include "includes/RMQRMM64.h"

using namespace std;
using namespace rmqrmm;

bool TRACE = false;			// true: print all details for console
bool TEST = true;			// true: apply exhaustive test
uint N_REP = 100000;		// number of REPETITIONS
bool RANDOM_TEST = true;	// true, random values (i,j) for test rmq->queryRMQ(i,j)
bool GROWDOWN = false;		// 1, increase pseudosorted arrays. 0 : decrease
bool LOAD_CREATE = true;	// 0/1 LOAD/CRATE
bool SEGMENT_EXP = false;
bool SAVElOAD_A = false;
bool isA = false;

int testRMQ(int argc, char *argv[]);
int testRMQLoad(int argc, char *argv[]);
double runExperimentTimeRMQ(RMQRMM64 *rmq, ulong n);
double runExperimentTimeRMQ_segmentSize(RMQRMM64 *rmq, ulong n, ulong seg, char fileName[400]);

int main(int argc, char *argv[]) {
	if(argc == 2){
		testRMQLoad(argc, argv);
	}else
		testRMQ(argc, argv);
}

int testRMQ(int argc, char *argv[]) {
	RMQRMM64 *rmq = NULL;
	uint pseudoSorted;
	char resultsFile[400], saveLoadFile[400], fileName[400], dirName[400];

	// different way to declare A array giving the different constructors of the RMQRMM64...
	//ulong *A;
	//short int *A;
	int *A;
	//long int *A;
	ulong n, i;
	double ct_rmq, qt_rmq;

	ct_rmq = 0.0;		// for construction time

	if(argc != 11){
		cout << "Missing parameters..." << endl;
		cout << "rmq_rmm_bp's usage: " << endl;
		cout << " " << argv[0] << " <array's length(n)> <load/create 0/1> <saveLoadFile> <repetitions> <pseudoSorted> <GROW/DOWN> <RandomWeight> <resultsFile> <dirName> <segment experiment>" << endl;
		return 1;
	}

	n = atol(argv[1]);
	LOAD_CREATE = atoi(argv[2]); // 0/1; 1: to create the structure from A ans store it in saveLoadFile. 0: to load the structure from saveLoadFile
	strcpy(saveLoadFile, argv[3]);
	N_REP = atol(argv[4]);
	pseudoSorted = atoi(argv[5]);
	GROWDOWN = atoi(argv[6]);
	ulong weight = atol(argv[7]);
	strcpy(resultsFile, argv[8]);
	strcpy(dirName, argv[9]);
	SEGMENT_EXP = atoi(argv[10]); // 0/1

	cout << "Parameters rmq_rmm_bp_3.0..."<< endl;
	cout << " Array's length(n) = " << n << endl;
	cout << " LOAD_CREATE = " << LOAD_CREATE << endl;
	cout << " saveLoadFile = " << saveLoadFile << endl;
	cout << " N_REP = " << N_REP << endl;
	cout << " pseudoSorted = " << pseudoSorted << endl;
	cout << " pseudoSorted GROW/DOWN = " << GROWDOWN << endl;
	cout << " weight = " << weight << endl;
	cout << " resultsFile = " << resultsFile << endl;
	cout << " dirName = " << dirName << endl;
	cout << " SEGMENT_EXP = " << SEGMENT_EXP << endl;

	// initialize random seed
	//srand(time(NULL));

	uint bpercell = ceilingLog64(weight+1, 2);
	cout << " bpercell = " << bpercell << endl;
	cout << " sizeof(ulong) " << sizeof(ulong) << endl;
	ulong num, cont, nA;
	nA = n*bpercell/(8*sizeof(ulong));
	if ((n*bpercell)%(8*sizeof(ulong)))
		nA++;

	if(LOAD_CREATE){ // create...
		//A = new ulong[nA];
		isA = true;
		//cout << " ** size of A array " << nA*sizeof(ulong) << " Bytes" << endl;
		//A = new short int[n];
		A = new int[n];
		//A = new long int[n];

		if(pseudoSorted){
			// IT WORKS ONLY FOR BIG INCREASING NUMBERS!
			if (GROWDOWN){
				cout << "Creating a Increasing sequence A[1.." << n << "]; A[i] is in [i-" << weight << ", i+ "<< weight << "]..." << endl;
				for(size_t i=0; i < n; i++)
					A[i] = (i-weight + rand()%(2*weight));
			}else{
				cout << "Creating a Decreasing sequence A[1.." << n << "]; A[i] is in [i-" << weight << ", i+ "<< weight << "]..." << endl;
				for(size_t i=0; i < n; i++)
					A[i] = n-i-weight + rand()%(2*weight);
			}
		}else{
			cout << "Creating the sequence A[1.." << n << "] with RANDOM numbers in the range [0.." << weight-1 << "]..." << endl;
			for(size_t i=0, cont=0; i < n; i++, cont+=bpercell){
				num = rand()%weight;
				A[i] = num;
				//setNum64(A, cont, bpercell, num);
			}
		}

		if(SAVElOAD_A){
			// Saving A[] for test in 'fileName'
			strcpy(fileName, "");
			strcpy(fileName, saveLoadFile);
			strcat(fileName, ".Seq");
			cout << "Saving A[] for test in " << fileName << endl;
			ofstream os (fileName, ios::binary);
			//nA = n*bpercell/8;
			//if((n*bpercell) % 8)	nA++;
			//os.write((const char*)A, nA*sizeof(ulong));
			//os.write((const char*)A, n*sizeof(short int));
			os.write((const char*)A, n*sizeof(int));
			//os.write((const char*)A, n*sizeof(long int));
			os.close();
		}

		if (TRACE){
			cout << "A[0.." << n-1 << "]" << endl;
			/*for (i=cont=0; i<n; i++, cont += bpercell){
				cout << getNum64(A, cont, bpercell) << " ";
				if ((i+1)%20==0)
					cout << endl;
			}
			cout << endl;*/
			for (i=0; i<n; i++){
				cout << A[i] << " ";
				if ((i+1)%20==0)
					cout << endl;
			}

		}

		cout << endl << "Creating rmq for A[]..." << endl;
		{
			ct_rmq = getTime_ms();
			//rmq = new RMQRMM64(A, bpercell, n, false);
			rmq = new RMQRMM64(A, n);
			ct_rmq = getTime_ms() - ct_rmq;
			//isA = false;
		}
		cout << "rmq creted in " << ct_rmq << " milliseconds" << endl;



		rmq->saveDS(saveLoadFile);
		rmq->~RMQRMM64();
		if(SAVElOAD_A) delete [] A;
		rmq->loadDS(saveLoadFile);
	}else	// load...
		rmq = new RMQRMM64(saveLoadFile);

	cout << "rmq size = " << rmq->getSize() << " Bytes = " << (float)rmq->getSize()*8.0/(float)n << " n" << endl;

	if (TEST){
		ulong n_min;
		cout << " Test RMQ..." << endl;

		if(SAVElOAD_A){
			// Loading A[] for test from 'fileName'
			//A = new ulong[nA];
			isA = true;
			//A = new short int[n];
			A = new int[n];
			//A = new long int[n];
			strcpy(fileName, "");
			strcpy(fileName, saveLoadFile);
			strcat(fileName, ".Seq");
			cout << "Loading A[] for test in " << fileName << endl;
			ifstream is(fileName, ios::binary);
			is.read((char*)A, nA*sizeof(ulong));
			//is.read((char*)A, n*sizeof(short int));
			//is.read((char*)A, n*sizeof(int));
			//is.read((char*)A, n*sizeof(long int));
			is.close();
		}

		if (!isA)
			cout << "It can not test... A[] array is NULL" << endl;
		else{
			if (TRACE){
				cout << "A[0.." << n-1 << "]" << endl;
				/*for (i=cont=0; i<n; i++, cont += bpercell){
					cout << getNum64(A, cont, bpercell) << " ";
					if ((i+1)%20==0)
						cout << endl;
				}
				cout << endl;*/
				for (i=0; i<n; i++){
					cout << A[i] << " ";
					if ((i+1)%20==0)
						cout << endl;
				}
			}

			ulong min, rmq_min, beg, end;
			beg=34;
			end=88;
			min = beg;
			for (i=beg+1; i<=end; i++){
				if (A[i] <= A[min])
					min = i;
			}
			rmq_min = rmq->queryRMQ(beg,end);
			if (rmq_min < beg || rmq_min > end){
				cout << "ERROR... rmq_min = " << rmq_min << " out of range [" << beg << " , " << end << " ]" << endl;
				exit(1);
			}else{
				if (A[rmq_min] != A[min]){
				//num = getNum64(A, rmq_min*bpercell, bpercell);
				//n_min = getNum64(A, min*bpercell, bpercell);
				//if (num != n_min){
					cout << "ERROR... rmq(" << beg << " , " << end << " ) = " << rmq_min << " != " << min << endl;
					exit(1);
				}
			}

			if (RANDOM_TEST){
				cout << "Random test O(N_REP)..." << endl;
				for (uint t=0; t<5000; t++){
					beg = (rand() % (n/2));
					end = n/2 + (rand() % (n/2)-1);
					if (end > beg+1){
						//n_min = getNum64(A, beg*bpercell, bpercell);
						for (min=beg, i=beg+1; i<=end; i++){
							if (A[i] <= A[min])
								min = i;

							//num = getNum64(A, i*bpercell, bpercell);
							//if(num < n_min){
								//min = i;
								//n_min = getNum64(A, min*bpercell, bpercell);
							//}
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
							//num = getNum64(A, rmq_min*bpercell, bpercell);
							//n_min = getNum64(A, min*bpercell, bpercell);
							//if (num != n_min){
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
						//num = getNum64(A, end*bpercell, bpercell);
						//n_min = getNum64(A, min*bpercell, bpercell);
						//if (num < n_min) min = end;

						if (A[end] < A[min]) min = end;

						rmq_min = rmq->queryRMQ(beg,end);
						if (rmq_min < beg || rmq_min > end){
							cout << "ERROR... rmq_min = " << rmq_min << " out of range [" << beg << " , " << end << " ]" << endl;
							exit(1);
						}else{
							if (A[rmq_min] != A[min]){
							//num = getNum64(A, rmq_min*bpercell, bpercell);
							//n_min = getNum64(A, min*bpercell, bpercell);
							//if (num != n_min){
								cout << "ERROR... (" << beg << " , " << end << " ) = " << rmq_min << " != " << min << endl;
								exit(1);
							}
						}
					}
				}
			}

			cout << "rmq test completed!" << endl;
		}
	}

	FILE *fp = fopen(resultsFile, "a+" );
	if(SEGMENT_EXP){
		for(ulong seg=10; seg <= 10000000; seg*=10){
			if(n > seg){
				strcpy(fileName, "");
				strcpy(fileName, dirName);
				qt_rmq = runExperimentTimeRMQ_segmentSize(rmq, n, seg, fileName);

				// [n] [o(n) extra size] [query time ms] [construction time ms] [segment] [weight] [BLK] [SB] [SS]
				fprintf(fp, "%lu %f %G %G %lu %lu %d %d\n", n, ((float)rmq->getSize()*8.0/(float)n)-2.0, qt_rmq, ct_rmq, seg, weight, BLK, SS);
			}
		}
	}else{
		qt_rmq = runExperimentTimeRMQ(rmq, n);
		// [n] [o(n) extra size] [query time ms] [construction time ms] [weight] [BLK] [SB] [SS]
		fprintf(fp, "%lu %f %G %G %lu %d %d\n", n, ((float)rmq->getSize()*8.0/(float)n)-2.0, qt_rmq, ct_rmq, weight, BLK, SS);
	}

	fclose(fp);

	rmq->~RMQRMM64();
	cout << "Program Finished OK" << endl;	return 0;
}


int testRMQLoad(int argc, char *argv[]) {
	cout << "Loading..." << endl;
	char* file = new char[100];
	strcpy(file, argv[1]);
	RMQRMM64 *rmq = new RMQRMM64(file);

	//cout << "rmq size = " << rmq->getSize() << " Bytes = " << (float)rmq->getSize()*8.0/((float)rmq->nP/2.0) << " n" << endl;
	runExperimentTimeRMQ(rmq, rmq->nP/2);
	rmq->~RMQRMM64();
	cout << "Program Finished OK" << endl;	return 0;
}


double runExperimentTimeRMQ_segmentSize(RMQRMM64 *rmq, ulong n, ulong seg, char fileName[400]){
	double t, avgTime;
	ulong i, j, k;
	ulong *QP = new ulong[N_REP];

	char *str = new char[100];
	strcpy(str, "");
	sprintf(str, "n%lus%lu.q", n, seg);
	strcat(fileName, str);

	cout << "_______________________________________________________________________________________________________" << endl;
	cout << "Start Experiment (random values for rmq(i,j) in fixed segment of length " << seg << endl;
	avgTime = 0.0;

	for (k=0; k<N_REP; k++){
		i = (rand() % (n-seg-1));
		j = i + seg;
		QP[k] = i;
		{
			t = getTime_ms();
			rmq->queryRMQ(i,j);
			t = getTime_ms() - t;
		}
		avgTime += t;
	}

	// Saving QP[]...
	cout << "Saving QP[] in " << fileName << endl;
	ofstream os (fileName, ios::binary);
	os.write((const char*)QP, N_REP*sizeof(unsigned long int));
	os.close();

	avgTime /= (double)N_REP;
	cout << "Average CPU time for item : " << avgTime*1000.0 << " Microseconds" << endl;
	cout << "_______________________________________________________________________________________________________" << endl;

	return avgTime*1000.0;
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

	return avgTime*1000.0;
}
