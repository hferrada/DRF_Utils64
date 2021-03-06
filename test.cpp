//============================================================================
// Name        : DRF_Utils64.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "includes/RMQRMM64.h"
#include "includes/ConfigFile.h"

using namespace std;
using namespace drf64;

bool TRACE = false;		// true: print all details for console
bool TEST = true;		// true: apply exhaustive test
uint N_REP = 10000;		// number of REPETITIONS
bool RANDOM = false;		// true, random vales (i,j) for test rmq->queryRMQ(i,j)

int testRMQRMM64(int argc, char *argv[]);
int testRMQRMM64Load(int argc, char *argv[]);
int testPropertyFile(int argc, char *argv[]);
void runExperimentTimeRMQ(RMQRMM64 *rmq, ulong n);

int main(int argc, char *argv[]) {

	if(argc == 2){
		testRMQRMM64Load(argc, argv);
	}else
		testRMQRMM64(argc, argv);
	//testPropertyFile(argc, argv);
}

int testRMQRMM64Load(int argc, char *argv[]) {
	cout << "Loading..." << endl;
	char* file = new char[100];
	strcpy(file, argv[1]);
	RMQRMM64 *rmq = new RMQRMM64(file);

	cout << "rmq size = " << rmq->getSize() << " Bytes = " << (float)rmq->getSize()*8.0/((float)rmq->nP/2.0) << " n" << endl;
	runExperimentTimeRMQ(rmq, rmq->nP/2);
	rmq->~RMQRMM64();
	cout << "Program Finished OK" << endl;	return 0;
}

int testRMQRMM64(int argc, char *argv[]) {
	long int *A;
	ulong n, i;

	if(argc != 3){
		cout << "Missing parameters..." << endl;
		cout << "RMQRMM64's usage: " << endl;
		cout << " " << argv[0] << " <array length(n)> <file name to save DT>" << endl;
		return 1;
	}

	n = atol(argv[1]);
	A = new long int[n];
	ulong weight = 1000;
	if (n > 100000)
		weight = 10000;
	if (n < 1000)
		weight = 100;

	cout << "Creating the sequence A[0.." << n-1 << "] with numbers in the range [0.." << weight-1 << "]..." << endl;
	for(i=0; i<n; i++)
		A[i] = (rand()%weight);

	if (TRACE){
		cout << "A[0.." << n-1 << "]" << endl;
		cout << A[0] << " ";
		for (i=1; i<n; i++){
			if (i%10 == 0)
				cout << " - ";
			cout << A[i] << " ";
		}
		cout << endl;
	}

	cout << "It is creating rmq for A[]..." << endl;
	RMQRMM64 *rmq = new RMQRMM64(A, n);
	char* file = new char[100];
	strcpy(file, argv[2]);
	rmq->saveDS(file);
	rmq->~RMQRMM64();
	rmq->loadDS(file);

	cout << "rmq size = " << rmq->getSize() << " Bytes = " << (float)rmq->getSize()*8.0/(float)n << " n" << endl;
	if (TEST){
		cout << "Test RMQ..." << endl;
		ulong min, rmq_min, beg, end;

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
	runExperimentTimeRMQ(rmq, n);
	rmq->~RMQRMM64();
	cout << "Program Finished OK" << endl;	return 0;
}

void runExperimentTimeRMQ(RMQRMM64 *rmq, ulong n){
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
}

int testPropertyFile(int argc, char *argv[]) {
	ConfigFile cf("config.txt");
	std::string foo;
	std::string water;
	double      four;

	foo   = (std::string)(cf.Value("section_1","foo"));
	water = (std::string)(cf.Value("section_2","water"));
	four  = (double)(cf.Value("section_2","four"));

	std::cout << foo   << std::endl;
	std::cout << water << std::endl;
	std::cout << four  << std::endl;

	return 0;
}
