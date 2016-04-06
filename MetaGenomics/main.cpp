/*
 * main.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider, Abhishek Biswas
 * Version: 2.0 (Alpha)
 */


#include "Common.h"
#include "Read.h"
#include "Dataset.h"
#include "HashTable.h"
#include "Edge.h"
#include "OverlapGraph.h"

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames,string & allFileName, UINT64 & minimumOverlapLength, bool & startFromUnitigGraph, UINT64 & maxThreads, UINT64 & writeGraphSize);

int main(int argc, char **argv)
{
	int numprocs, myid;
	double start, end;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
	start = MPI_Wtime();

	UINT64 minimumOverlapLength;
	vector<string> pairedEndFileNames, singleEndFileNames;
	string allFileName;
	bool startFromUnitigGraph = false;
	UINT64 maxThreads = DEF_THREAD_COUNT;
	UINT64 writeGraphSize = MAX_PAR_GRAPH_SIZE;
	parseArguments(argc, argv, pairedEndFileNames, singleEndFileNames, allFileName, minimumOverlapLength, startFromUnitigGraph, maxThreads, writeGraphSize);
	Dataset *dataSet = new Dataset(pairedEndFileNames, singleEndFileNames, allFileName, minimumOverlapLength);
	HashTable *hashTable=new HashTable(numprocs);
	hashTable->insertDataset(dataSet, minimumOverlapLength,numprocs, myid);
	//OverlapGraph *overlapGraph;
	//overlapGraph=new OverlapGraph(hashTable,maxThreads,writeGraphSize,allFileName); //hashTable deleted by this function after building the graph also writes graph
	delete hashTable;	// Do not need the hash table any more.
	//delete dataSet;
	//delete overlapGraph;

	MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
	end = MPI_Wtime();

	if (myid == 0) { /* use time on master node */
	    printf("Runtime for %d processes = %f\n", numprocs, end-start);
	}
	MPI_Finalize();
}

/**********************************************************************************************************************
	Parse the input arguments
**********************************************************************************************************************/

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames, string & allFileName, UINT64 & minimumOverlapLength, bool & startFromUnitigGraph, UINT64 & maxthreads, UINT64 & writeGraphSize)
{
	allFileName = "";
	minimumOverlapLength = 0;
	startFromUnitigGraph = false;
	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for(int i = 0; i < argc; i++)
	{
		cout << argv[i] << ' ';
	}
	cout << endl;
	while(argc--)
			argumentsList.push_back(*argv++);

	if(argumentsList.size() == 1)
	{
		cerr << endl << "Usage: MetaGenomics [OPTION]...[PRARAM]..." << endl;
		cerr << "  -pe\tnumber of files and paired-end file names" <<endl; 			// Paired-end file name in fasta/fastq format. mate pairs should be one after another in the file.
		cerr << "  -se\tnumber of files and single-end file names" <<endl; 			// Single-end file name in fasta/fastq format.
		cerr << "  -f\tAll file name prefix" <<endl; 			// all output file with have this name with different extensions.
		cerr << "  -l\tminimum overlap length" << endl; 	// Minimum overlap length for two reads to overlap in the overlap graph.
		cerr << "  -t\tmaximum threads used" << endl; 	// Maximum OMP threads used
		cerr << "  -w\tgraph write frequency per" << endl; 	// Maximum size of sub-graph before its written to disk
		cerr << "  -s\tstart from unitig graph" << endl; 	// -s means that the program will build the graph. Otherwise it will load the graph from the unitig graph file.
			exit(0);
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		if(argumentsList[i] == "-pe")
		{
			UINT64 numberOfPairedEndDatasets = atoi(argumentsList[++i].c_str());
			for(UINT64 j = 0; j < numberOfPairedEndDatasets; j++)
			{
				pairedEndFileNames.push_back(argumentsList[++i]);
			}
		}
		else if(argumentsList[i] == "-se")
		{
			UINT64 numberOfSingleEndDatasets = atoi(argumentsList[++i].c_str());
			for(UINT64 j = 0; j < numberOfSingleEndDatasets; j++)
			{
				singleEndFileNames.push_back(argumentsList[++i]);
			}
		}
		else if (argumentsList[i] == "-f")
			allFileName = argumentsList[++i];
		else if (argumentsList[i] == "-l")
			minimumOverlapLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-s")
			startFromUnitigGraph = true;
		else if (argumentsList[i] == "-t")
			maxthreads = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-w")
			writeGraphSize = atoi(argumentsList[++i].c_str());
		else
		{
			cerr << endl << "Usage: MetaGenomics [OPTION]...[PRARAM]..." << endl;
			cerr << "  -pe\tnumber of files and paired-end file names" <<endl; 			// Paired-end file name in fasta/fastq format. mate pairs should be one after another in the file.
			cerr << "  -se\tnumber of files and single-end file names" <<endl; 			// Single-end file name in fasta/fastq format.
			cerr << "  -f\tAll file name prefix" <<endl; 			// all output file with have this name with different extensions.
			cerr << "  -l\tminimum overlap length" << endl; 	// Minimum overlap length for two reads to overlap in the overlap graph.
			cerr << "  -w\tgraph write frequency per" << endl; 	// Maximum size of sub-graph before it's written to disk
			cerr << "  -t\tmaximum threads used" << endl; 	// Maximum OMP threads used
			cerr << "  -s\tstart from unitig graph" << endl; 	// -s means that the program will build the graph. Otherwise it will load the graph from the unitig graph file.
			if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
				exit(0);
			else
			{
				cerr << "Unknown option: " << argumentsList[i] << endl << endl;
				exit(1);
			}
		}
	}
}
