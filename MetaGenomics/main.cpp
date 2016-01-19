/*
 * main.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 * Version: 0.1 (Alpha)
 */


#include "Common.h"
#include "Read.h"
#include "Dataset.h"
#include "HashTable.h"
#include "Edge.h"
#include "OverlapGraph.h"





void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames,string & allFileName, UINT64 & minimumOverlapLength, bool & startFromUnitigGraph);

int main(int argc, char **argv)
{
	CLOCKSTART;
	UINT64 minimumOverlapLength;
	vector<string> pairedEndFileNames, singleEndFileNames;
	UINT64 counter;
	UINT64 iteration = 0;
	string allFileName;
	bool startFromUnitigGraph = false;
	parseArguments(argc, argv, pairedEndFileNames, singleEndFileNames, allFileName, minimumOverlapLength, startFromUnitigGraph);
	Dataset *dataSet = new Dataset(pairedEndFileNames, singleEndFileNames, minimumOverlapLength);
	OverlapGraph *overlapGraph;

	if(startFromUnitigGraph) // Read the graph from unitig file
	{
		overlapGraph = new OverlapGraph();
		overlapGraph->setDataset(dataSet);
		overlapGraph->readGraphFromFile(allFileName+".unitig");
		overlapGraph->sortEdges();
	}
	else // Build the graph from the scratch
	{
		HashTable *hashTable=new HashTable();
		hashTable->insertDataset(dataSet, minimumOverlapLength);
		overlapGraph=new OverlapGraph(hashTable); //hashTable deleted by this function after building the graph
		dataSet->saveReads(allFileName+"_sortedReads.fasta");
		overlapGraph->sortEdges();
		overlapGraph->saveGraphToFile(allFileName+".unitig");
	}

	overlapGraph->calculateFlow(allFileName+"_flow.input", allFileName+"_flow.output");
	cout << "nodes: " << overlapGraph->getNumberOfNodes() << " edges: " << overlapGraph->getNumberOfEdges() << endl;
	overlapGraph->printGraph(allFileName+"graph1.gdl", allFileName+"contigs1.fasta");

	overlapGraph->removeAllSimpleEdgesWithoutFlow();

	overlapGraph->calculateMeanAndSdOfInsertSize();



	do
	{
		// Mate pair paths are used to simplify the graph in this step
		cout << endl;
		cout << "===============================================================================================================================================" <<endl;
		cout << "FIRST LOOP ITERATION " << ++iteration << endl;
		cout << "===============================================================================================================================================" <<endl;
		overlapGraph->simplifyGraph();
		counter = overlapGraph->findSupportByMatepairsAndMerge();
	} while (counter > 0 && iteration < loopLimit); // To avoid infinite loops

	overlapGraph->printGraph(allFileName+"graph2.gdl", allFileName+"contigs2.fasta");

	iteration = 0;
	do
	{
		// Scaffolder
		cout << endl;
		cout << "===============================================================================================================================================" <<endl;
		cout << "SECOND LOOP ITERATION " << ++iteration << endl;
		cout << "===============================================================================================================================================" <<endl;
		overlapGraph->simplifyGraph();
		counter = overlapGraph->scaffolder();

	} while (counter > 0 && iteration < loopLimit);// To avoid infinite loops

	overlapGraph->printGraph(allFileName+"graph3.gdl", allFileName+"contigs3.fasta");

	iteration = 0;
	do
	{
		// Coverage depth information is used to resolve ambiguity
		cout << endl;
		cout << "===============================================================================================================================================" <<endl;
		cout << "THIRD LOOP ITERATION " << ++iteration << endl;
		cout << "===============================================================================================================================================" <<endl;
		overlapGraph->simplifyGraph();
		counter = overlapGraph->resolveNodes();

	} while (counter > 0 && iteration < loopLimit);// To avoid infinite loops

	overlapGraph->printGraph(allFileName+"graph4.gdl", allFileName+"contigs4.fasta");

	delete dataSet;
	delete overlapGraph;
	CLOCKSTOP;
}



/**********************************************************************************************************************
	Parse the input arguments
**********************************************************************************************************************/

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames, string & allFileName, UINT64 & minimumOverlapLength, bool & startFromUnitigGraph)
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
		else
		{
			cerr << endl << "Usage: MetaGenomics [OPTION]...[PRARAM]..." << endl;
			cerr << "  -pe\tnumber of files and paired-end file names" <<endl; 			// Paired-end file name in fasta/fastq format. mate pairs should be one after another in the file.
			cerr << "  -se\tnumber of files and single-end file names" <<endl; 			// Single-end file name in fasta/fastq format.
			cerr << "  -f\tAll file name prefix" <<endl; 			// all output file with have this name with different extensions.
			cerr << "  -l\tminimum overlap length" << endl; 	// Minimum overlap length for two reads to overlap in the overlap graph.
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
