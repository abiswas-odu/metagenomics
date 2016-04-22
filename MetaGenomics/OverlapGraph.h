/*
 * OverlapGraph.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_
#include "Common.h"
#include "Dataset.h"
#include "HashTable.h"
#include "Edge.h"

/**********************************************************************************************************************
	Class to store the overlap graph.
**********************************************************************************************************************/

#define MIN_MARKED 128				//No. of reads to be marked before a send communication is initiated...

enum nodeType {
	UNEXPLORED = 0, // Current node u is not explored yet. Meaning that there is no edge (u,v) in the graph.
	EXPLORED = 1, //  Current node u is explored. Meaning that all edges (u,v) are inserted in the dataset.
	EXPLORED_AND_TRANSITIVE_EDGES_MARKED = 2, // Meaning that all transitive edges (u,v) of current node u is marked and its neighbors transitive edges are also marked. Now it is safe to remove any transitive edge (u,v) from node u.
	EXPLORED_AND_TRANSITIVE_EDGES_REMOVED = 3, // Meaning that all transitive edges (u,v) of current node u has been removed and its neighbors transitive edges are also marked. Now it is safe to remove any transitive edge (u,v) from node u.
	EXPLORED_AND_TRANSITIVE_EDGES_WRITTEN = 4
};

enum markType{
	VACANT = 0,
	INPLAY = 1,
	ELIMINATED = 2
};

class OverlapGraph
{
	private:
		Dataset * dataSet; 											// Pointer to the dataset containing all the reads.
		HashTable * hashTable;										// Pointer to the hash table.
		vector<UINT64> meanOfInsertSizes; 							// Mean of insert sizes.
		vector<UINT64> sdOfInsertSizes; 							// Standard deviation of insert sizes.
		UINT64 estimatedGenomeSize;									// Estimated genome size. Works for isolated genome. Will not work for Metagenomics.
		int myProcID;													// Id of the MPI process
		UINT64 numberOfNodes;										// Number of nodes in the overlap graph.
		UINT64 numberOfEdges;										// Number of edges in the overlap graph.
		UINT64 parallelThreadPoolSize;								//No. of OMP threads to spawn
		UINT64 writeParGraphSize;									//No. of vertices to mark before writing graph to memory
		INT64 longestMeanOfInsertSize; // CP: the longest mean insert size out of all datasets: max(meanOfInsertSizes[i])
		UINT8 mergedEdgeOrientation(Edge *edge1, Edge *edge2);		// Orientation of the edge when two edges are merged.
		UINT8 twinEdgeOrientation(UINT8 orientation);				// Orientation of the reverse edge.
		bool mergeList(Edge *edge1, Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * ListOrientations);
		bool findPathBetweenMatepairs(Read * read1, Read * read2, UINT8 orient, UINT8 datasetNumbe, vector <Edge *> &copyOfPath, vector <UINT64> &copyOfFlags);
		UINT64 exploreGraph(Edge* firstEdge, Edge * lastEdge, UINT64 distanceOnFirstEdge, UINT64 distanceOnLastEdge, UINT64 datasetNumber, UINT64 level, vector <Edge *> &firstPath, vector <UINT64> &flags);
	public:
		bool flowComputed;											// Flag to check wheather the flow is computed or not.
		OverlapGraph(void);											// Default constructor.
		OverlapGraph(HashTable *ht,UINT64 maxThreads,UINT64 maxParGraph, string fnamePrefix,int myid, int MPINodeBlockSize, int numprocs);								// Another constructor.
		~OverlapGraph();											// Destructor.
		bool markTransitiveEdges(UINT64 readNumber, map<UINT64, vector<Edge*> * > *parGraph); // Mark transitive edges of a read.
		bool buildOverlapGraphFromHashTable(HashTable *ht, string fnamePrefix, int MPINodeBlockSize, int numprocs);			// Build the overlap graph using hashtable.
		bool insertEdge(Edge * edge, map<UINT64, vector<Edge*> * > *parGraph); 								// Insert an edge in the partial overlap graph.
		bool insertEdge(Read *read1, Read *read2, UINT64 r1Len, UINT64 r2Len,  UINT8 orient, UINT16 overlapOffset, map<UINT64, vector<Edge*> * > *parGraph); // Insert an edge in the overlap graph.

		bool checkOverlapForContainedRead(string read1, string read2, UINT64 orient, UINT64 start);
		bool checkOverlap(string read1, string read2, UINT64 orient, UINT64 start);

		//bool checkOverlap(Read *read1, Read *read2, UINT64 orient, UINT64 start); // Check overlap between two reads after a match is found using the hash table.
		//bool checkOverlapForContainedRead(Read *read1, Read *read2, UINT64 orient, UINT64 start);
		bool insertAllEdgesOfRead(UINT64 readNumber, map<UINT64,nodeType> * exploredReads, map<UINT64, vector<Edge*> * > *parGraph);	// Insert into the overlap graph all edges of a read.
		bool removeTransitiveEdges(UINT64 readNumber, map<UINT64, vector<Edge*> * > *parGraph);				// Remove all transitive edges from the overlap graph incident to a given read.
		bool saveParGraphToFile(string fileName, map<UINT64,nodeType> * exploredReads, map<UINT64, vector<Edge*> * > *parGraph);   //Save partial graph to file and reduce memory footprint
		void markContainedReads(string fnamePrefix, int numprocs);									// Find superReads for each read and mark them as contained read.
};



#endif /* OVERLAPGRAPH_H_ */
