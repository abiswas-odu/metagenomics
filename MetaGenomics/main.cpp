/*
 * main.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider, Abhishek Biswas
 * Version: 3.o
 */


#include "Common.h"
#include "Read.h"
#include "Dataset.h"
#include "HashTable.h"
#include "Edge.h"
#include "OverlapGraph.h"

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames,string & allFileName, UINT64 & minimumOverlapLength, bool & startFromUnitigGraph, UINT64 & maxThreads, UINT64 & writeGraphSize, UINT64 &maxMemSizeGB);

int main(int argc, char **argv)
{
	CLOCKSTART;
	UINT64 minimumOverlapLength;
	vector<string> pairedEndFileNames, singleEndFileNames;
	string allFileName;
	bool startFromUnitigGraph = false;
	UINT64 maxThreads = DEF_THREAD_COUNT;
	UINT64 writeGraphSize = MAX_PAR_GRAPH_SIZE;
	UINT64 maxMemSizeGB = getMaxMemory();
	cout<<"Max available memory: "<<maxMemSizeGB<< " GB"<<endl;
	parseArguments(argc, argv, pairedEndFileNames, singleEndFileNames, allFileName, minimumOverlapLength, startFromUnitigGraph, maxThreads, writeGraphSize,maxMemSizeGB);
	cout<<"Max usable memory: "<<maxMemSizeGB<< " GB"<<endl;
	Dataset *dataSet = new Dataset(pairedEndFileNames, singleEndFileNames, allFileName, minimumOverlapLength);
	OverlapGraph *overlapGraph;
	HashTable *hashTable=new HashTable();
	hashTable->insertDataset(dataSet, minimumOverlapLength,maxThreads);
	overlapGraph=new OverlapGraph(hashTable,maxThreads,writeGraphSize,allFileName); //hashTable deleted by this function after building the graph also writes graph
	delete dataSet;
	delete overlapGraph;
	CLOCKSTOP;
}

/**********************************************************************************************************************
	Parse the input arguments
**********************************************************************************************************************/

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames, string & allFileName, UINT64 & minimumOverlapLength, bool & startFromUnitigGraph, UINT64 & maxthreads, UINT64 & writeGraphSize, UINT64 &maxMemSizeGB)
{
	allFileName = "";
	minimumOverlapLength = 0;
	startFromUnitigGraph = false;
	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for(UINT64 i = 0; i < argc; i++)
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
		cerr << "  -m\tmaximum memory usage allowed" << endl; 	// Maximum memory to be used before written to disk
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
		else if (argumentsList[i] == "-m")
			maxMemSizeGB = atoi(argumentsList[++i].c_str());
		else
		{
			cerr << endl << "Usage: MetaGenomics [OPTION]...[PRARAM]..." << endl;
			cerr << "  -pe\tnumber of files and paired-end file names" <<endl; 			// Paired-end file name in fasta/fastq format. mate pairs should be one after another in the file.
			cerr << "  -se\tnumber of files and single-end file names" <<endl; 			// Single-end file name in fasta/fastq format.
			cerr << "  -f\tAll file name prefix" <<endl; 			// all output file with have this name with different extensions.
			cerr << "  -l\tminimum overlap length" << endl; 	// Minimum overlap length for two reads to overlap in the overlap graph.
			cerr << "  -w\tgraph write frequency per" << endl; 	// Maximum size of sub-graph before it's written to disk
			cerr << "  -t\tmaximum threads used" << endl; 	// Maximum OMP threads used
			cerr << "  -m\tmaximum memory usage allowed" << endl; 	// Maximum memory to be used before written to disk
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
