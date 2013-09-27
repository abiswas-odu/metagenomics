/*
 * OverlapGraph.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#include "Common.h"
#include "OverlapGraph.h"
#include "CS2/cs2.h"


/**********************************************************************************************************************
	Check if two edges match.
	e1(u,v) and e2(v,w). At node v, one of the edges should be an incoming edge and the other should be an outgoing
	edge to match.
**********************************************************************************************************************/

bool matchEdgeType(const Edge *edge1, const Edge *edge2)
{
	if     ( (edge1->getOrientation() == 1 || edge1->getOrientation() == 3) && (edge2->getOrientation() == 2 || edge2->getOrientation() == 3) ) // *-----> and >------*
		return true;
	else if( (edge1->getOrientation() == 0 || edge1->getOrientation() == 2) && (edge2->getOrientation() == 0 || edge2->getOrientation() == 1) ) // *------< and <------*
		return true;
	return false;
}


/**********************************************************************************************************************
	Function to compare two edges. Used for sorting.
**********************************************************************************************************************/
bool compareEdgeID (const Edge *edge1, const Edge* edge2)
{
	return (edge1->getDestinationRead()->getReadNumber() < edge2->getDestinationRead()->getReadNumber());
}

/**********************************************************************************************************************
	Function to compare two edges. Used for sorting.
**********************************************************************************************************************/
// CP: Sort by overlap offset. What for?
// BH: In the transitive edge removal step, we need the edges to be sorted according to their overlap offset.
bool compareEdges (const Edge *edge1, const Edge* edge2)
{
	return (edge1->getOverlapOffset() < edge2->getOverlapOffset());
}

// CP: if we want to use a different way to decide to resolve a cross by coverage depth, we just need to change this function, right?
// BH: Yes you only need to change this function to match the coverage depth.
bool isOverlappintInterval(INT64 mean1, INT64 sd1, INT64 mean2, INT64 sd2)
{
	int start1 = mean1 - 2*sd1;
	int end1 = mean1 + 2*sd2;
	int start2 = mean2 - 2*sd2;
	int end2 = mean2 + 2 *sd2;
	return (start1 >= start2 && start1 <= end2) || (end1 >= start2 && end1 <= end2) || (start2 >= start1 && start2 <= end1) || (end2 >= start1 && end2 <= end1);
}

/**********************************************************************************************************************
	Default Constructor
**********************************************************************************************************************/
OverlapGraph::OverlapGraph(void)
{
	// Initialize the variables.
	estimatedGenomeSize = 0;
	numberOfNodes = 0;
	numberOfEdges = 0;
	flowComputed = false;
}



/**********************************************************************************************************************
	Another Constructor. Build the overlap grpah using the hash table.
**********************************************************************************************************************/
/*
OverlapGraph::OverlapGraph(HashTable *ht)
{
	// Initialize the variables.
	estimatedGenomeSize = 0;
	numberOfNodes = 0;
	numberOfEdges = 0;
	flowComputed = false;
	buildOverlapGraphFromHashTable(ht);
}
*/



/**********************************************************************************************************************
	Default destructor.
**********************************************************************************************************************/
OverlapGraph::~OverlapGraph()
{
	// Free the memory used by the overlap graph.
	for(UINT64 i = 0; i < graph->size(); i++)
	{
		for(UINT64 j = 0; j< graph->at(i)->size(); j++)
		{
			delete graph->at(i)->at(j);
		}
		delete graph->at(i);
	}
	delete graph;
}


/**********************************************************************************************************************
	Build the overlap graph from hash table
**********************************************************************************************************************/
bool OverlapGraph::buildOverlapGraphFromHashTable(HashTable *ht)
{
	CLOCKSTART;
	estimatedGenomeSize = 0;
	numberOfNodes = 0;
	numberOfEdges = 0;
	flowComputed = false;
	hashTable = ht;
	dataSet = ht->getDataset();
	UINT64 counter = 0;
	vector<nodeType> *exploredReads = new vector<nodeType>;
	exploredReads->reserve(dataSet->getNumberOfUniqueReads()+1);

	// This is a queue used to process the reads. The reads are processed in a certain order so that we can remove the transitive edges after inserting all the edges of a read.
	// Once all edges of a read is inserted, then we add all the edges of neighbors and then neighbors' neighbors.
	vector<UINT64> * queue = new vector<UINT64>;
	queue->reserve(dataSet->getNumberOfUniqueReads()+1);

	//This vector is used for discovering transitive edges.
	vector<markType> *markedNodes = new vector<markType>;
	markedNodes->reserve(dataSet->getNumberOfUniqueReads()+1);

	// CP: Initialize graph containing every unique read represented by a node.
	graph = new vector< vector<Edge *> * >;
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);



	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // Initialization
	{
		vector<Edge *> *newList = new vector<Edge *>;
		graph->push_back(newList);
		exploredReads->push_back(UNEXPLORED);
		queue->push_back(0);
		markedNodes->push_back(VACANT);
	}
	// CP: assign contained reads to super-reads (populate superReadID in the Read class), but not remove them
	markContainedReads();

	// CP: why reading the mate pair information here, not when constructing Dataset initially?
	// BH: We read the matepair information here again. Because some of the reads are contained and will be assigned to their super read ID in the mate-pair.
	this->dataSet->readMatePairsFromFile();

	// CP: find the overlap between reads and remove transitive edges along the way
	for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(exploredReads->at(i) == UNEXPLORED)
		{
			UINT64 start = 0, end = 0; 											// Initialize queue start and end.
			queue->at(end) = i;													// CP: use the value of end [i.e. queue->at(0) = i], THEN increment end.
			end ++;
			while(start < end) 													// This loop will explore all connected component starting from read i.
			{
				counter++;
				UINT64 read1 = queue->at(start++);
				if(exploredReads->at(read1) == UNEXPLORED)
				{
					insertAllEdgesOfRead(read1, exploredReads);					// Explore current node.
					exploredReads->at(read1) = EXPLORED;
				}
				if(graph->at(read1)->size() != 0) 								// Read has some edges (required only for the first read when a new queue starts.
				{
					if(exploredReads->at(read1) == EXPLORED) 					// Explore unexplored neighbors first.
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++ )
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead()->getReadNumber();
							if(exploredReads->at(read2) == UNEXPLORED) 			// Not explored.
							{
								queue->at(end++) = read2; 						// Put in the queue.
								insertAllEdgesOfRead(read2, exploredReads);
								exploredReads->at(read2) = EXPLORED;
							}
						}
						markTransitiveEdges(read1, markedNodes); // Mark transitive edges
						exploredReads->at(read1) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
					}
					if(exploredReads->at(read1) == EXPLORED_AND_TRANSITIVE_EDGES_MARKED)
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++) 				// Then explore all neighbour's neighbors
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead()->getReadNumber();
							if(exploredReads->at(read2) == EXPLORED)
							{
								for(UINT64 index2 = 0; index2 < graph->at(read2)->size(); index2++) 		// Explore all neighbors neighbors
								{
									UINT64 read3 = graph->at(read2)->at(index2)->getDestinationRead()->getReadNumber();
									if(exploredReads->at(read3) == UNEXPLORED) 				// Not explored
									{
										queue->at(end++) = read3; 					// Put in the queue
										insertAllEdgesOfRead(read3, exploredReads);
										exploredReads->at(read3) = EXPLORED;
									}
								}
								markTransitiveEdges(read2, markedNodes); // Mark transitive edge
								exploredReads->at(read2) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
							}
						}
						removeTransitiveEdges(read1); // Remove the transitive edges
					}
				}
				if(counter%100000==0)	// Show the progress.
					cout<<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges/2 << endl;
			}
		}
	}
	cout<<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges/2 << endl;

	delete exploredReads;
	delete queue;
	delete markedNodes;
	//delete hashTable;	// Do not need the hash table any more.
	do
	{
		 counter = contractCompositePaths();
		 counter += removeDeadEndNodes();
	} while (counter > 0);

	// CP: nodes that don't have any edge will be ignored in all operations
	// BH: If a node u does not have any edge in the graph then the list of edges of node u will be empty in the graph.

	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
	This function check if a read contains other small reads. If a read is contained in more than one super read
	then it is assigned to the longest such super read.
**********************************************************************************************************************/
void OverlapGraph::markContainedReads(void)
{
	CLOCKSTART;
	if(dataSet->longestReadLength == dataSet->shortestReadLength) // If all reads are of same length, then no need to do look for contained reads.
	{
		cout << "All reads are of same length. No contained reads." << endl;
		CLOCKSTOP;
		return;
	}
	UINT64 counter = 0;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read
	{
		Read *read1 = dataSet->getReadFromID(i); // Get the read
		string readString = read1->getStringForward(); // Get the forward of the read
		string subString;
		for(UINT64 j = 1; j < read1->getReadLength() - hashTable->getHashStringLength(); j++) // fGr each substring of read1 of length getHashStringLength
		{
			subString = readString.substr(j,hashTable->getHashStringLength()); // Get the substring from read1
			vector<UINT64> * listOfReads=hashTable->getListOfReads(subString); // Search the substring in the hash table
			if(!listOfReads->empty()) // If other reads contain the substring as prefix or suffix
			{
				for(UINT64 k = 0; k < listOfReads->size(); k++) // For each read in the list.
				{
					UINT64 data = listOfReads->at(k); // We used bit operation in the hash table to store read ID and orientation
					Read *read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
																						// Most significant 2 bits store the orientation.
																						// Orientation 0 means prefix of forward of the read
																						// Orientation 1 means suffix of forward of the read
																						// Orientation 2 means prefix of reverse of the read
																						// Orientation 3 means prefix of reverse of the read

					if(readString.length() > read2->getStringForward().length() && checkOverlapForContainedRead(read1,read2,(data >> 62),j)) // read1 need to be longer than read2 in order to contain read2
																																			 // Check if the remaining of the strings also match
					{
						if(read2->superReadID == 0) // This is the rist super read found. we store the ID of the super read.
						{
							read2->superReadID = i;
							counter ++;
						}
						else // Already found some previous super read. Now we have another super read.
						{
							if(readString.length() > dataSet->getReadFromID(read2->superReadID)->getReadLength()) // This super read is longer than the previous super read. Update the super read ID.
								read2->superReadID = i;
						}
					}
				}
			}
		}
		if(i%1000000 == 0)
			cout << setw(10) << counter << " contained reads in " << setw(10) << i << " super reads." << endl;
	}

	// Get some statistics
	UINT64 containedReads = 0, nonContainedReads = 0;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		Read *rr = dataSet->getReadFromID(i);
		if(rr->superReadID == 0) // Count the number of reads that are not contained by some other reads.
			nonContainedReads++;
		else					// Count the number of reads that are contained by some other read.
			containedReads++;
	}
	cout<< endl << setw(10) << nonContainedReads << " Non-contained reads. (Keep as is)" << endl;
	cout<< setw(10) << containedReads << " contained reads. (Need to change their mate-pair information)" << endl;
	CLOCKSTOP;
}


/**********************************************************************************************************************
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depends on the orient).
	orient 0 means prefix of forward of the read2
	orient 1 means suffix of forward of the read2
	orient 2 means prefix of reverse of the read2
	orient 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read2 is contained in read1.
**********************************************************************************************************************/
bool OverlapGraph::checkOverlapForContainedRead(const Read *read1, const Read *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getStringForward(); // Get the forward of read1
	UINT64 hashStringLength = hashTable->getHashStringLength(), lengthRemaining1, lengthRemaining2;
	string string2 = (orient == 0 || orient== 1) ? read2->getStringForward() : read2->getStringReverse(); // Get the string in read2 based on the orientation.
	if(orient == 0 || orient == 2)
									// orient 0
									//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
									//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
									//				OR
									// orient 2
									//	 >---*****MMMMMMMMMMMMMMM*******------> read1
									//		      MMMMMMMMMMMMMMM*******<	    Reverse complement of read2
	{
		lengthRemaining1 = string1.length() - start - hashStringLength; 	// This is the remaining of read1
		lengthRemaining2 = string2.length() - hashStringLength; 	// This is the remaining of read2
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return string1.substr(start + hashStringLength, lengthRemaining2) == string2.substr(hashStringLength, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	else							// orient 1
									//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
									//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
									//				OR
									// orient 3
									//	 >---*****MMMMMMMMMMMMMMM-------------> read1
									//		<*****MMMMMMMMMMMMMMM				Reverse Complement of Read2
	{
		lengthRemaining1 = start;
		lengthRemaining2 = string2.length() - hashStringLength;
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return string1.substr(start - lengthRemaining2, lengthRemaining2) == string2.substr(0, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	return false;

}



/**********************************************************************************************************************
	Checks if two read overlaps.
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	Orientation 0 means prefix of forward of the read2
	Orientation 1 means suffix of forward of the read2
	Orientation 2 means prefix of reverse of the read2
	Orientation 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read1 and read2 overlap.
**********************************************************************************************************************/
bool OverlapGraph::checkOverlap(const Read *read1, const Read *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getStringForward(); // Get the forward string of read1
	UINT64 hashStringLength = hashTable->getHashStringLength();
	string string2 = (orient == 0 || orient== 1) ? read2->getStringForward() : read2->getStringReverse(); // Get the string from read2 according to orient.
	if(orient == 0 || orient == 2)		// orient 0
										//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
										//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
										//				OR
										// orient 2
										//	 >---*****MMMMMMMMMMMMMMM*************> 			read1
										//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2
	{
		if(string1.length()- start - hashStringLength >= string2.length() - hashStringLength) // The overlap must continue till the end.
			return false;
		return string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)) == string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)); // If the remaining strings match.
	}
	else								// orient 1
										//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
										//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
										//				OR
										// orient 3
										//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
										//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2
	{
		if(string2.length()-hashStringLength < start)
			return false;
		return string1.substr(0, start) == string2.substr(string2.length()-hashStringLength-start, start); // If the remaining strings match.
	}
}



/**********************************************************************************************************************
	Insert an edge in the overlap graph.
**********************************************************************************************************************/
bool OverlapGraph::insertEdge(Edge * edge)
{

	UINT64 ID = edge->getSourceRead()->getReadNumber(); // This is the source read.
	if(graph->at(ID)->empty()) 							// If there is no edge incident to the node
		numberOfNodes++;								// Then a new node is inserted in the graph. Number of nodes increased.
	graph->at(ID)->push_back(edge);						// Insert the edge in the list of edges of ID
		numberOfEdges++;								// Increase the number of edges.
	updateReadLocations(edge);							// If the current edge contains some reads, then we need to update their location information.
	return true;
}



/**********************************************************************************************************************
	Insert an edge in the graph.
**********************************************************************************************************************/
bool OverlapGraph::insertEdge(Read *read1, Read *read2, UINT8 orient, UINT16 overlapOffset)
{
	Edge * edge1 = new Edge(read1,read2,orient,overlapOffset);								// Create a new edge in the graph to insert.
	UINT16 overlapOffsetReverse = read2->getReadLength() + overlapOffset - read1->getReadLength();	// Set the overlap offset accordingly for the reverse edge. Note that read lengths are different.
																						// If read lengths are the same. Then the reverse edge has the same overlap offset.
	Edge * edge2 = new Edge(read2,read1,twinEdgeOrientation(orient),overlapOffsetReverse);		// Create a new edge for the reverses string.

	edge1->setReverseEdge(edge2);		// Set the reverse edge pointer.
	edge2->setReverseEdge(edge1);		// Set the reverse edge pinter.
	insertEdge(edge1);					// Insert the edge in the overlap graph.
	insertEdge(edge2);					// Insert the edge in the overlap graph.
	return true;
}


/**********************************************************************************************************************
	This function prints the overlap graph in overlap_graph->gdl file. The graph can be viewed by
	aisee (free software available at http://www.aisee.com/)
	It also stores the contigs in a file.


graph: {
layoutalgorithm :forcedir
fdmax:704
tempmax:254
tempmin:0
temptreshold:3
tempscheme:3
tempfactor:1.08
randomfactor:100
gravity:0.0
repulsion:161
attraction:43
ignore_singles:yes
node.fontname:"helvB10"
edge.fontname:"helvB10"
node.shape:box
node.width:80
node.height:20
node.borderwidth:1
node.bordercolor:31
node: { title:"43" label: "43" }	// node, Title and label are both node ID 43 (and read number)
node: { title:"65" label: "65" }
............................................
............................................
// edges from source node 43 to destination node 32217, thickness of 3 means composite edge, thickness of 1 for simple edge
// edge type of backarrowstyle:solid arrowstyle:solid color: green is >----------------<
// edge type of arrowstyle:solid color: red is <----------------<
// edge type of arrowstyle: none color: blue  is <------------------->
// (1,0x,206,30) means (Flow, coverageDepth, OverlapOffset, numberOfReads)

edge: { source:"43" target:"32217" thickness: 3 backarrowstyle:solid arrowstyle:solid color: green label: "(1,0x,206,30)" }
edge: { source:"65" target:"38076" thickness: 3 arrowstyle:solid color: red label: "(0,0x,75,11)" }
edge: { source:"280" target:"47580" thickness: 3 arrowstyle: none color: blue label: "(0,0x,123,11)" }
}

**********************************************************************************************************************/
bool OverlapGraph::printGraph(string graphFileName, string contigFileName)
{
	CLOCKSTART;
	vector <Edge *> contigEdges;
	UINT64 thickness, highestDegree = 0, highestDegreeNode;
	ofstream graphFilePointer, contigFilePointer;
	graphFilePointer.open(graphFileName.c_str());
	contigFilePointer.open(contigFileName.c_str());
	if(graphFilePointer == NULL)
		MYEXIT("Unable to open file: "+graphFileName);
	if(contigFilePointer == NULL)
		MYEXIT("Unable to open file: "+contigFileName);
	graphFilePointer << "graph: {" << endl <<  "layoutalgorithm :forcedir" << endl <<  "fdmax:704" << endl <<  "tempmax:254" << endl <<  "tempmin:0" << endl <<  "temptreshold:3" << endl <<  "tempscheme:3" << endl <<  "tempfactor:1.08" << endl <<  "randomfactor:100" << endl <<  "gravity:0.0" << endl <<  "repulsion:161" << endl <<  "attraction:43" << endl <<  "ignore_singles:yes" << endl <<  "node.fontname:\"helvB10\"" << endl << "edge.fontname:\"helvB10\"" << endl <<  "node.shape:box" << endl <<  "node.width:80" << endl <<  "node.height:20" << endl <<  "node.borderwidth:1" << endl <<  "node.bordercolor:31" << endl;
	for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
		if(!graph->at(i)->empty())
			graphFilePointer << "node: { title:\""<<i<<"\" label: \"" << i << "\" }" << endl;


	for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
	{
			if(!graph->at(i)->empty())
			{
				if(graph->at(i)->size() > highestDegree)
				{
					highestDegree = graph->at(i)->size();
					highestDegreeNode = i;

				}
				for(UINT64 j=0; j < graph->at(i)->size(); j++)
				{
					Edge * e = graph->at(i)->at(j);
					UINT64 source = e->getSourceRead()->getReadNumber(), destination = e->getDestinationRead()->getReadNumber();
					if(source < destination || (source == destination && e < e->getReverseEdge()) )
					{
						contigEdges.push_back(e); // List of contigs.
						thickness = e->getListOfReads()->empty() ? 1: 3;
						if(e->getOrientation() == 0)
							graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " arrowstyle: none backarrowstyle: solid color: red label: \"(" << e->flow << "," << e->coverageDepth << "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
						else if(e->getOrientation() == 1)
							graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " backarrowstyle:solid arrowstyle:solid color: green label: \"(" << e->flow << "," << e->coverageDepth <<  "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
						else if(e->getOrientation() == 2)
							graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " arrowstyle: none color: blue label: \"(" << e->flow << "," << e->coverageDepth << "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
						else if(e->getOrientation() == 3)
							graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " arrowstyle:solid color: red label: \"(" << e->flow << "," << e->coverageDepth <<  "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
					}
				}
			}
	}
	graphFilePointer << "}";
	cout << "Aisee graph written." << endl;
	sort(contigEdges.begin(),contigEdges.end(),compareEdges); // Sort the contigs by their length.
	reverse(contigEdges.begin(), contigEdges.end());
	UINT64 sum = 0;
	for(UINT64 i = 0; i < contigEdges.size(); i++) // Store the contigs in a file.
	{
		string s = getStringInEdge(contigEdges.at(i)); // get the string in the edge.
		contigFilePointer << ">contig_"<< i+1 << " Flow: " << setw(10) << contigEdges.at(i)->flow <<  " Edge  (" << setw(10) << contigEdges.at(i)->getSourceRead()->getReadNumber() << ", " << setw(10) << contigEdges.at(i)->getDestinationRead()->getReadNumber() << ") String Length: " << setw(10) << s.length() << " Coverage: " << setw(10) << contigEdges.at(i)->coverageDepth << endl;
		sum += s.length();
		UINT64 start=0;
		do
		{
			contigFilePointer << s.substr(start,100) << endl;  // save 100 BP in each line.
			start+=100;
		} while (start < s.length());
	}

	// Print some statistics about the graph.
	cout<< "Total contig length: " << sum <<" BP" << endl;
	cout<< "Number of Nodes in the graph: " << getNumberOfNodes() << endl;
	cout<< "Number of Edges in the graph: " << getNumberOfEdges()/2 << endl;

	UINT64 simEdges = 0, comEdges = 0, inEdges = 0, outEdges = 0;
	for(UINT64 i=0; i < graph->at(highestDegreeNode)->size(); i++)
	{
		if(graph->at(highestDegreeNode)->at(i)->getListOfReads()->empty())
			simEdges++;
		else
			comEdges++;
		if(graph->at(highestDegreeNode)->at(i)->getOrientation() == 0 || graph->at(highestDegreeNode)->at(i)->getOrientation() == 1)
			inEdges++;
		else
			outEdges++;
	}
	// Print some more statistics on the node with highest degree.
	cout<< "Highest Degree Read " << highestDegreeNode << " has " << highestDegree << " neighbors." << endl;
	cout << "In Edges: " << inEdges << " Out Edges: " << outEdges << " Simple Edges: " << simEdges << " Composite Edges: " << comEdges << endl;
	cout<< "String: " << dataSet->getReadFromID(highestDegreeNode)->getStringForward() << endl;

	graphFilePointer.close();
	contigFilePointer.close();
	CLOCKSTOP;
	return true;
}





/**********************************************************************************************************************
	Insert all edges of a read in the overlap graph
**********************************************************************************************************************/
bool OverlapGraph::insertAllEdgesOfRead(UINT64 readNumber, vector<nodeType> * exploredReads)
{
	Read *read1 = dataSet->getReadFromID(readNumber); 	// Get the current read read1.
	string readString = read1->getStringForward(); 		// Get the forward string of read1.
	string subString;
	for(UINT64 j = 1; j < read1->getReadLength()-hashTable->getHashStringLength(); j++) // For each proper substring of length getHashStringLength of read1
	{
		subString = readString.substr(j,hashTable->getHashStringLength());  // Get the proper substring s of read1.
		vector<UINT64> * listOfReads=hashTable->getListOfReads(subString); // Search the string in the hash table.
		if(!listOfReads->empty()) // If there are some reads that contain s as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT16 overlapOffset;
				UINT8 orientation;
				Read *read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				if(exploredReads->at(read2->getReadNumber())!= UNEXPLORED) 			// No need to discover the same edge again. All edges of read2 is already inserted in the graph.
						continue;
				if(read1->superReadID == 0 && read2->superReadID == 0 && checkOverlap(read1,read2,(data >> 62),j)) // Both read need to be non contained.
				{
					switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
						case 0: orientation = 3; overlapOffset = read1->getReadLength() - j; break; 				// 3 = r1>------->r2
						case 1: orientation = 0; overlapOffset = hashTable->getHashStringLength() + j; break; 		// 0 = r1<-------<r2
						case 2: orientation = 2; overlapOffset = read1->getReadLength() - j; break; 				// 2 = r1>-------<r2
						case 3: orientation = 1; overlapOffset = hashTable->getHashStringLength() + j; break; 		// 1 = r2<------->r2
					}
					insertEdge(read1,read2,orientation,read1->getStringForward().length()-overlapOffset); 			// Insert the edge in the graph.
				}
			}
		}
	}
	if(graph->at(readNumber)->size() != 0)
		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdges); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}




/**********************************************************************************************************************
	Mark all the transitive edges of a read.
	For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
**********************************************************************************************************************/
bool OverlapGraph::markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes)
{

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++) // Mark all the neighbours of the current read as INPLAY
		markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber()) = INPLAY; // Inplay

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++) // Traverse through the list of edges according to their overlap offset.
	{
		UINT64 read2 = graph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber(); // For each neighbor
		if(markedNodes->at(read2) == INPLAY) 										// If the neighbor is marked as INPLAY
		{
			for(UINT64 j = 0; j < graph->at(read2)->size(); j++)
			{
				UINT64 read3 = graph->at(read2)->at(j)->getDestinationRead()->getReadNumber(); // Get the neighbors neighbors
				if(markedNodes->at(read3) == INPLAY)
				{

					UINT8 type1 = graph->at(readNumber)->at(i)->getOrientation();
					UINT8 type2 = graph->at(read2)->at(j)->getOrientation();
					if((type1 == 0 ||  type1 == 2) && (type2==0 || type2==1)) 	// Check edge orientation
						markedNodes->at(read3) = ELIMINATED; 					// Mark as ELIMINATED
					else if((type1==1||type1==3) && (type2==2 || type2==3)) 	// Check edge orientation
						markedNodes->at(read3) = ELIMINATED; 					// Mark as ELIMINATED
				}
			}
		}
	}
	for(UINT64 i = 0;i < graph->at(readNumber)->size(); i++)
	{
		if(markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber()) == ELIMINATED) // Current read to a node marked as ELIMINATED
		{
			graph->at(readNumber)->at(i)->transitiveRemovalFlag = true; 					// Mark this edge as transitive edge. Will remove this edge later.
			graph->at(readNumber)->at(i)->getReverseEdge()->transitiveRemovalFlag = true;	// Mark also the reverse edge. Will remove this edge later.
		}
	}

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++)
		markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber()) = VACANT; 	// Change back all the variables modified in this function to VACANT

	markedNodes->at(readNumber) = VACANT; 		// Mark as vacant.
	return true;
}



/**********************************************************************************************************************
	Remove all transitive edges of a given read.
	For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
**********************************************************************************************************************/
bool OverlapGraph::removeTransitiveEdges(UINT64 readNumber)
{
	for(UINT64 index = 0; index < graph->at(readNumber)->size(); index++)  		// Go through the list of edges of the current read.
	{
		if(graph->at(readNumber)->at(index)->transitiveRemovalFlag == true)		// This edge is marked as transitive. We will first remove the reverse edge.
		{
			Edge *twinEdge = graph->at(readNumber)->at(index)->getReverseEdge();
			UINT64 ID = twinEdge->getSourceRead()->getReadNumber();
			for(UINT64 index1 = 0; index1 < graph->at(ID)->size(); index1++) 	// Get the reverse edge first
			{
				if(graph->at(ID)->at(index1) == twinEdge)
				{
					delete twinEdge;
					graph->at(ID)->at(index1) = graph->at(ID)->at(graph->at(ID)->size()-1); // Move the transitive edges at the back of the list and remove.
					graph->at(ID)->pop_back();
					if(graph->at(ID)->empty())
						numberOfNodes--;
					numberOfEdges--;
					break;
				}
			}
		}
	}
	UINT64 j=0;
	for(UINT64 index=0; index < graph->at(readNumber)->size(); index++) // Then we will remove all the transitive edges of the current read.
	{
		if(graph->at(readNumber)->at(index)->transitiveRemovalFlag == false)		// We move all the non-transitive edges at the beginning of the list
			graph->at(readNumber)->at(j++) = graph->at(readNumber)->at(index);
		else		// Free the transitive edge
		{
			numberOfEdges--;
			delete graph->at(readNumber)->at(index);
		}
	}
	graph->at(readNumber)->resize(j);
	if(graph->at(readNumber)->empty())
		numberOfNodes--;
	return true;
}


/**********************************************************************************************************************
	Contract composite paths in the overlap graph.
	u*-------->v>---------*w  => u*--------------------*w
	u*--------<v<---------*w  => u*--------------------*w
**********************************************************************************************************************/
UINT64 OverlapGraph::contractCompositePaths(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	for(UINT64 index = 1 ; index < graph->size(); index++)
	{
		if(graph->at(index)->size() == 2) // Check if the node has only two edges.
		{
			Edge *edge1 = graph->at(index)->at(0);  // First edge.
			Edge *edge2 = graph->at(index)->at(1);  // Second Edge.
			if(flowComputed == true || !isEdgePresent(edge1->getDestinationRead()->getReadNumber(), edge2->getDestinationRead()->getReadNumber()))
																						// Before flow is computed we do not insert multiple edges between the same endpoints.
				// CP: Before flow is computed (flowComputed == false), why checking if there is an edge between the two destination reads? why don't you need to do this after flow is calculated?
				// BH: Before the flow is computed we do not want to insert multiple edges between the same nodes. This condition is required for CS2 minimum Cost flow algorithm
				// CP: Why the destination reads, not the source reads??
				// BH: We do not want to remove loops (a,a).o you remove
				// CP: using the example above, don't you need to check if u and w are different reads or not?
				// BH: We do not need to check if u and w are different read or not.
			{
				if( matchEdgeType(edge1->getReverseEdge(), edge2) && edge1->getSourceRead() != edge1->getDestinationRead()) // One incoming edge and one outgoing edge.
				{
					mergeEdges(edge1->getReverseEdge(),edge2);	// Merge the edges.
					counter++;									// Counter how many edges merged.
				}
			}
		}

	}
	cout << setw(10) << counter << " composite Edges merged." << endl;
	CLOCKSTOP;
	return counter;
}




/**********************************************************************************************************************
	Merge two edges in the overlap graph.
	CP: the flow of the new edge is the minimum flow of the two old edges and the flow of the new edge is deducted from those of the old edges
	CP: Remove the old edges that doesn't have any flow left, but keep the old edges that still have flow left.
**********************************************************************************************************************/
bool OverlapGraph::mergeEdges(Edge *edge1, Edge *edge2)
{
	Edge *edgeForward = new Edge(); // New forward edge.
	Read *read1 = edge1->getSourceRead(), *read2 = edge2->getDestinationRead();
	Edge *edgeReverse = new Edge(); // New reverse edge.

	UINT8 orientationForward = mergedEdgeOrientation(edge1,edge2);			// Orientation of the forward edge.
	UINT8 orientationReverse = twinEdgeOrientation(orientationForward);		// Orientation of the reverse edge.

	vector<UINT64> * listReadsForward = new vector<UINT64>;					// List of reads in the forward edge.
	vector<UINT16> * listOverlapOffsetsForward= new vector<UINT16>;				// List of Overlaps in the forward edge.
	vector<UINT8> * listOrientationsForward = new vector<UINT8>;			// List of orientations in the forward edge.

	mergeList(edge1, edge2, listReadsForward, listOverlapOffsetsForward, listOrientationsForward); // Merge the lists from the two edges.

	edgeForward->makeEdge(read1,read2,orientationForward, edge1->getOverlapOffset() + edge2->getOverlapOffset(), listReadsForward, listOverlapOffsetsForward, listOrientationsForward); // Make the forward edge

	vector<UINT64> * listReadsReverse = new vector<UINT64>;					// List of reads in the reverse edge.
	vector<UINT16> * listOverlapOffsetsReverse= new vector<UINT16>;				// List of overlaps in the reverse edge.
	vector<UINT8> * listOrientationsReverse = new vector<UINT8>;			// List of orientations in the reverse edge.

	mergeList(edge2->getReverseEdge(),edge1->getReverseEdge(), listReadsReverse, listOverlapOffsetsReverse,listOrientationsReverse);	// Merge the lists from the two reverse edges.

	edgeReverse->makeEdge(read2, read1, orientationReverse, edge2->getReverseEdge()->getOverlapOffset() + edge1->getReverseEdge()->getOverlapOffset(), listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse); // Make the reverse edge

	edgeForward->setReverseEdge(edgeReverse);		// Update the reverse edge pointer.
	edgeReverse->setReverseEdge(edgeForward);		// Update the reverse edge pointer.


	UINT16 flow = min(edge1->flow,edge2->flow);		// We take the minimum of the flow in the new edge.

	edgeForward->flow = flow;						// Modify the flow in the new forward edge.
	edgeReverse->flow = flow;						// Modify the flow in the new reverse edge.


	insertEdge(edgeForward);						// Insert the new forward edge in the graph.
	insertEdge(edgeReverse);						// Insert the new reverse edge in the graph.

	edge1->flow = edge1->flow - flow;				// Remove the used flow from edge1.
	edge1->getReverseEdge()->flow = edge1->flow;	// Remove the used flow from the reverse of edge1.

	edge2->flow = edge2->flow - flow;				// Remove the used flow from edge2.
	edge2->getReverseEdge()->flow = edge2->flow;	// Remove the used flow from teh reverse of edge2.

	if(edge1->flow == 0 || flow == 0)				// If no flow left in edge1
		removeEdge(edge1);							// edge1 is deleted from the graph.
	if(edge2->flow == 0 || flow == 0)				// If now flow left in edge2
		removeEdge(edge2);							// edge 2 is deleted from the graph.

	return true;

}


/**********************************************************************************************************************
	Merge the list of reads, list of overlap offsets and list of orientations of two edges.
**********************************************************************************************************************/

bool OverlapGraph::mergeList(const Edge *edge1, const Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations)
{
	UINT64 sum = 0;
	for(UINT64 i = 0; i < edge1->getListOfOrientations()->size(); i++) 			// Take the list from edge1.
	{
		listReads->push_back(edge1->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge1->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge1->getListOfOrientations()->at(i));
		sum += edge1->getListOfOverlapOffsets()->at(i);
	}
	listReads->push_back(edge1->getDestinationRead()->getReadNumber()); 		// Insert the common node of the two edges

	listOverlapOffsets->push_back(edge1->getOverlapOffset() - sum);					// Get the overlap offset.

	if(edge1->getOrientation() == 1 || edge1->getOrientation() == 3)			// Orientation of the common node. Depends on the orientation of the edges.
		listOrientations->push_back(1);
	else
		listOrientations->push_back(0);
	for(UINT64 i = 0; i < edge2->getListOfOrientations()->size(); i++)			// take the list from edge2.
	{
		listReads->push_back(edge2->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge2->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge2->getListOfOrientations()->at(i));
	}
	return true;
}




/**********************************************************************************************************************
	Orientation of the merged edge.
	e1(u,v) e2(v,w)
	Orientation e1 0 = u<-------<v		Orientation e2 0 = v<-------<w
	Orientation e1 1 = u<------->v		Orientation e2 1 = v<------->w
	Orientation e1 2 = u>-------<v		Orientation e2 2 = v>-------<w
	Orientation e1 3 = u>------->v		Orientation e2 3 = v>------->w

	0 + 0 = u<-------<w
	0 + 1 = u<------->w
	...................

**********************************************************************************************************************/
UINT8 OverlapGraph::mergedEdgeOrientation(const Edge *edge1, const Edge *edge2)
{
	UINT8 or1 = edge1->getOrientation(), or2 = edge2->getOrientation(),returnValue;
	if(or1 == 0 && or2 == 0)
		returnValue = 0;
	else if(or1 == 0 && or2 == 1)
		returnValue = 1;
	else if(or1 == 1 && or2 == 2)
		returnValue = 0;
	else if(or1 == 1 && or2 == 3)
		returnValue = 1;
	else if(or1 == 2 && or2 == 0)
		returnValue = 2;
	else if(or1 == 2 && or2 == 1)
		returnValue = 3;
	else if(or1 == 3 && or2 == 2)
			returnValue = 2;
	else if(or1 == 3 && or2 == 3)
			returnValue = 3;
	else
	{
		cout<<(int)or1<<" "<<(int)or2<<endl;
		MYEXIT("Unable to merge.")
	}
	return returnValue;
}





/**********************************************************************************************************************
	Orientation of a reverse edge;
	Twin edge of Orientation 0 = <-------< is Orientation 3 = >------->
	Twin edge of Orientation 1 = <-------> is Orientation 1 = <------->
	Twin edge of Orientation 2 = >-------< is Orientation 2 = >-------<
	Twin edge of Orientation 3 = >-------> is Orientation 0 = <-------<
**********************************************************************************************************************/
UINT8 OverlapGraph::twinEdgeOrientation(UINT8 orientation)
{
	UINT8 returnValue;
	if(orientation == 0)
		returnValue = 3;
	else if(orientation == 1)
		returnValue = 1;
	else if(orientation == 2)
		returnValue = 2;
	else if(orientation == 3)
		returnValue = 0;
	else
		MYEXIT("Unsupported edge orientation.")
	return returnValue;
}




/**********************************************************************************************************************
	remove an edge from the overlap graph.
**********************************************************************************************************************/
bool OverlapGraph::removeEdge(Edge *edge)
{
	removeReadLocations(edge);								// If the current edge contains some reads. We have to update their location formation.
	removeReadLocations(edge->getReverseEdge());			// Same for the reverse edge.
	Edge *twinEdge = edge->getReverseEdge();
	UINT64 ID1 = edge->getSourceRead()->getReadNumber(), ID2 = edge->getDestinationRead()->getReadNumber();  // Get the source and destation read IDs.
	for(UINT64 i = 0; i< graph->at(ID2)->size(); i++) // Delete the twin edge first.
	{
		if(graph->at(ID2)->at(i) == twinEdge)
		{
			delete graph->at(ID2)->at(i);
			graph->at(ID2)->at(i) = graph->at(ID2)->at(graph->at(ID2)->size()-1);
			graph->at(ID2)->pop_back();
			if(graph->at(ID2)->empty())
				numberOfNodes--;
			numberOfEdges--;
			break;
		}
	}
	for(UINT64 i = 0; i< graph->at(ID1)->size(); i++) // Delete the edge then.
	{
		if(graph->at(ID1)->at(i) == edge)
		{
			delete graph->at(ID1)->at(i);
			graph->at(ID1)->at(i) = graph->at(ID1)->at(graph->at(ID1)->size()-1);
			graph->at(ID1)->pop_back();
			if(graph->at(ID1)->empty())
				numberOfNodes--;
			numberOfEdges--;
			break;
		}
	}
	return true;
}



/**********************************************************************************************************************
	Remove an all simple edge in the overlap graph that does not have any flow.
	CP: return the number of edges removed
**********************************************************************************************************************/
UINT64 OverlapGraph::removeAllSimpleEdgesWithoutFlow()
{
	CLOCKSTART;
	vector <Edge *> listOfEdges;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read.
	{
		if(!graph->at(i)->empty())	// If the read has some edges.
		{
			for(UINT64 j=0; j < graph->at(i)->size(); j++) // For each edge
			{
				Edge * edge = graph->at(i)->at(j);
			if(edge->getSourceRead()->getReadNumber() < edge->getDestinationRead()->getReadNumber() && edge->getListOfReads()->empty() && edge->flow == 0) // The edge is simple edge with no flow.
				{
					listOfEdges.push_back(edge); // Put in the list of edges to be removed.
				}
			}
		}
	}
	for(UINT64 i = 0 ; i < listOfEdges.size(); i++)
		removeEdge(listOfEdges.at(i));		// remove the edges from the list.
	cout<<"Edges removed: " << listOfEdges.size() << endl;
	CLOCKSTOP;
	return listOfEdges.size();
}


/**********************************************************************************************************************
	Remove an all edge in the overlap graph that does not have any flow.
	CP: return the number of edges removed
**********************************************************************************************************************/
UINT64 OverlapGraph::removeAllEdgesWithoutFlow()
{
	CLOCKSTART;
	vector <Edge *> listOfEdges;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read.
	{
		if(!graph->at(i)->empty())	// If the read has some edges.
		{
			for(UINT64 j=0; j < graph->at(i)->size(); j++) // For each edge
			{
				Edge * edge = graph->at(i)->at(j);
			if(edge->getSourceRead()->getReadNumber() < edge->getDestinationRead()->getReadNumber() && edge->flow == 0) // The edge is simple edge with no flow.
				{
					listOfEdges.push_back(edge); // Put in the list of edges to be removed.
				}
			}
		}
	}
	for(UINT64 i = 0 ; i < listOfEdges.size(); i++)
		removeEdge(listOfEdges.at(i));		// remove the edges from the list.
	cout<<"Edges removed: " << listOfEdges.size() << endl;
	CLOCKSTOP;
	return listOfEdges.size();
}


/**********************************************************************************************************************
	Remove nodes with all simple edges and all same arrow type

	A node is all incomming edges or all outgoing edges is called a dead end node. While traversing the graph if we
	enter such node, there is no way we can go out. So we remove such nodes from the graph. To remove the node all
	of its edges must be simple edges or very short edges (less than 10 read in it).
**********************************************************************************************************************/
UINT64 OverlapGraph::removeDeadEndNodes(void)
{
	CLOCKSTART;
	vector <UINT64> listOfNodes; // for saving nodes that should be deleted
	UINT64 edgesRemoved = 0;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read.
	{
		if(!graph->at(i)->empty())	// If the read has some edges.
		{
			// CP: comment on what's flag, inEdge and outEdge. It's hard to figure them out by reading the functions below
			// BH: flag is used to check if we need to remove the current node.
			// Initailly flag is set to 0. If the current node has an edge with more than 10 reads in it, then the flag will be 1 and the node will not be removed.
			UINT64 flag = 0, inEdge = 0 , outEdge = 0;
			for(UINT64 j=0; j < graph->at(i)->size(); j++) // For each edge
			{
				Edge * edge = graph->at(i)->at(j);
				// Break case:
				// 1. if composite edge is more than deadEndLength reads
				// 2. if the edge is loop for the current node
				// Then flag=1 and exit the loop
				// If the current node has an edge with more than 10 reads in it then we do not mark it as dead end.
				if(edge->getListOfReads()->size() > deadEndLength || edge->getSourceRead()->getReadNumber() == edge->getDestinationRead()->getReadNumber())
				{
					flag = 1;
					break;
				}

				if(edge->getOrientation() == 0 || edge->getOrientation() == 1)
					inEdge++;			// the number of incoming edges to this node
				else
					outEdge++;			// the number of outgoing edges to this node
			}
			// if all edges are incoming or outgoing AND has less than 10 reads, this node is marked to be removed?
			// This node does not have any composite edge with more than 10 reads in it.
			if(flag == 0) // If not break case
			{
				if( (inEdge > 0 && outEdge == 0) || (inEdge == 0 && outEdge > 0)) // only one type of simple edges
				{
					listOfNodes.push_back(i);
				}
			}
		}
	}

	//  in below actually to remove the marked deadend nodes and all their edges.
	vector <Edge *> listOfEdges;
	for(UINT64 i = 0 ; i < listOfNodes.size(); i++)
	{
		listOfEdges.clear();
		if(!graph->at(listOfNodes.at(i))->empty())	// If the read has some edges.
		{
			edgesRemoved += graph->at(listOfNodes.at(i))->size();
			for(UINT64 j=0; j < graph->at(listOfNodes.at(i))->size(); j++) // For each edge
			{
				listOfEdges.push_back(graph->at(listOfNodes.at(i))->at(j));
			}
			for(UINT64 j = 0; j< listOfEdges.size(); j++)
				removeEdge(listOfEdges.at(j));							// Remove all the edges of the current node.
		}
	}
	cout<< "Dead-end nodes removed: " << listOfNodes.size() << endl;
	cout<<  "Total Edges removed: " << edgesRemoved << endl;
	CLOCKSTOP;
	return listOfNodes.size();
}



/**********************************************************************************************************************
	Estimate the genome size. The estimated genome size is not used anywere in the program.
**********************************************************************************************************************/
bool OverlapGraph::estimateGenomeSize(void)
{
	CLOCKSTART;
	UINT64 delta, deltaSum, freqSum, currentGenomeSize = 0, previousGenomeSize = 0, freq, counter = 0;
	do
	{
		counter++; 						// To avoid infinite loop.
		deltaSum = 0;
		freqSum = 0;
		float a_statistic;
		for(UINT64 i = 1; i < graph->size(); i++)
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{
				Edge * currentEdge = graph->at(i)->at(j);
				if(currentEdge->getSourceRead()->getReadNumber() < currentEdge->getDestinationRead()->getReadNumber())
				{
					delta = currentEdge->getOverlapOffset();
					freq = 0;
					for(UINT64 k = 0; k < currentEdge->getListOfReads()->size(); k++)
						freq +=  dataSet->getReadFromID( currentEdge->getListOfReads()->at(k) )->getFrequency();
					if(previousGenomeSize != 0) 					// If not the first estimation
					{
						a_statistic=(float)((float)delta*(float)((float)dataSet->getNumberOfReads() / (float)previousGenomeSize) - (float)((float)freq * log((float)2)));
						if(a_statistic >= aStatisticsThreshold && delta >= minDelta)
						{
							deltaSum += delta;
							freqSum += freq;
						}
					}
					else if(currentEdge->getOverlapOffset() > 500)  // Only take edges longer than 500 for the first estimation
					{
						deltaSum += delta;
						freqSum += freq;
					}
				}
			}
		}
		previousGenomeSize = currentGenomeSize;
		currentGenomeSize = (int)((float)dataSet->getNumberOfReads()/(float)freqSum*(float)deltaSum);
		cout<<"Current estimated genome size: " << currentGenomeSize << endl;
	}while (currentGenomeSize != previousGenomeSize && counter < 10);
	estimatedGenomeSize=currentGenomeSize;
	cout<<"Final estimated genome size: " << getEstimatedGenomeSize() << endl;
	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
	Update location of reads for the new edge.
	The new edge many contain some reads. We need to update the location information of all such read.
**********************************************************************************************************************/
bool OverlapGraph::updateReadLocations(Edge *edge)
{
	UINT64 distance = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) 				// For each read in the edge
	{
		distance += edge->getListOfOverlapOffsets()->at(i);					// Distance of the read in the edge.
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));	// Get the read object.
		if(edge->getListOfOrientations()->at(i) == 1)						// Orientation 1 means that the forward string of the read is contained in this edge.
		{
			read->getListOfEdgesForward()->push_back(edge);					// Insert the current edge in the list of forward edges.
			read->getLocationOnEdgeForward()->push_back(distance);			// Also insert the distance within the edge.
			read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
			read->getLocationOnEdgeForward()->resize(read->getLocationOnEdgeForward()->size());
		}
		else																// Orientation 0 means that the reverser string of the read is contained in this edge.
		{
			read->getListOfEdgesReverse()->push_back(edge);					// Insert the current edge in the list of reverser edges.
			read->getLocationOnEdgeReverse()->push_back(distance);			// Also insert the distance withing the edge.
			read->getListOfEdgesReverse()->resize(read->getListOfEdgesReverse()->size());
			read->getLocationOnEdgeReverse()->resize(read->getLocationOnEdgeReverse()->size());
		}
	}
	return true;
}




/**********************************************************************************************************************
	Remove read mapping information of the current edge. This function is called when an edge is removed.
**********************************************************************************************************************/
bool OverlapGraph::removeReadLocations(Edge *edge)
{
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) 						// For each read in this edge.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i)); 		// Get the read object.
		for(UINT64 j = 0; j < read->getListOfEdgesForward()->size(); j++)  			// Check the list of edges that contain forward of this read.
		{
			if(read->getListOfEdgesForward()->at(j) == edge)						// Remove the current edge from the list.
			{
				read->getListOfEdgesForward()->at(j) = read->getListOfEdgesForward()->at(read->getListOfEdgesForward()->size()-1);
				read->getLocationOnEdgeForward()->at(j) = read->getLocationOnEdgeForward()->at(read->getLocationOnEdgeForward()->size()-1);

				read->getListOfEdgesForward()->pop_back();
				read->getLocationOnEdgeForward()->pop_back();

				read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
				read->getLocationOnEdgeForward()->resize(read->getLocationOnEdgeForward()->size());
			}
		}

		for(UINT64 j = 0; j < read->getListOfEdgesReverse()->size(); j++) 			// Check the list of edges that contain reverse of this read.
		{
			if(read->getListOfEdgesReverse()->at(j) == edge) 						// Remove the current edge from the list.
			{
				read->getListOfEdgesReverse()->at(j) = read->getListOfEdgesReverse()->at(read->getListOfEdgesReverse()->size()-1);
				read->getLocationOnEdgeReverse()->at(j) = read->getLocationOnEdgeReverse()->at(read->getLocationOnEdgeReverse()->size()-1);

				read->getListOfEdgesReverse()->pop_back();
				read->getLocationOnEdgeReverse()->pop_back();

				read->getListOfEdgesReverse()->resize(read->getListOfEdgesReverse()->size());
				read->getLocationOnEdgeReverse()->resize(read->getLocationOnEdgeReverse()->size());
			}
		}
	}
	return true;
}





/**********************************************************************************************************************
	Calculate Mean and Standard Deviation of insert size of each dataset.
**********************************************************************************************************************/
bool OverlapGraph::calculateMeanAndSdOfInsertSize(void)
{
	CLOCKSTART;
	longestMeanOfInsertSize = 0;

	vector<INT64> *insertSizes = new vector<INT64>;
	vector<Edge *> listOfEdgesRead1, listOfEdgesRead2;
	vector<INT64> locationOnEdgeRead1, locationOnEdgeRead2;

	for(UINT64 d = 0; d < dataSet->numberOfPairedDatasets; d++)	// For each dataset.
	{

		cout << "Calculating mean and SD of dataset: " << d << endl;
		insertSizes->clear();

		for(UINT64 i = 1; i < graph->size(); i++)	// for each read in the dataset
		{
			Read *read1 = dataSet->getReadFromID(i), *read2;	// Get a read read1 in the graph.

			for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++) 	// For each mate-pair read2
			{
				if(read1->getMatePairList()->at(j).datasetNumber == d)		// if it is in the current dataset.
				{
					listOfEdgesRead1.clear();  locationOnEdgeRead1.clear();
					listOfEdgesRead1.insert(listOfEdgesRead1.end(), read1->getListOfEdgesForward()->begin(), read1->getListOfEdgesForward()->end());				// All the edges that contain forward string of read1
					listOfEdgesRead1.insert(listOfEdgesRead1.end(), read1->getListOfEdgesReverse()->begin(), read1->getListOfEdgesReverse()->end());				// All the edges that contain reverse string of read1

					locationOnEdgeRead1.insert(locationOnEdgeRead1.end(), read1->getLocationOnEdgeForward()->begin(), read1->getLocationOnEdgeForward()->end());  	// Location on the corresponding edge.
					locationOnEdgeRead1.insert(locationOnEdgeRead1.end(), read1->getLocationOnEdgeReverse()->begin(), read1->getLocationOnEdgeReverse()->end());	// Location on the corresponding edge.

					read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);	// Get the read2 object.
					listOfEdgesRead2.clear();  locationOnEdgeRead2.clear();

					listOfEdgesRead2.insert(listOfEdgesRead2.end(), read2->getListOfEdgesForward()->begin(), read2->getListOfEdgesForward()->end());				// All the edges that contain forward string of read2
					listOfEdgesRead2.insert(listOfEdgesRead2.end(), read2->getListOfEdgesReverse()->begin(), read2->getListOfEdgesReverse()->end());				// All the edges that contain reverse string of read1

					locationOnEdgeRead2.insert(locationOnEdgeRead2.end(), read2->getLocationOnEdgeForward()->begin(), read2->getLocationOnEdgeForward()->end());	// Location on the corresponding edge.
					locationOnEdgeRead2.insert(locationOnEdgeRead2.end(), read2->getLocationOnEdgeReverse()->begin(), read2->getLocationOnEdgeReverse()->end());	// Location on the corresponding edge.

					for(UINT64 k = 0; k < listOfEdgesRead1.size(); k++)				// For each edge containing read1
					{
						for(UINT64 l = 0; l < listOfEdgesRead2.size(); l++)			// For each edge containing read2
						{
							if(listOfEdgesRead1.at(k) == listOfEdgesRead2.at(l) && locationOnEdgeRead1.at(k) > locationOnEdgeRead2.at(l))		// Both reads are on the same edge
							{
									if(locationOnEdgeRead1.at(k) - locationOnEdgeRead2.at(l) < 1000)					// Distance between the two edges is less than 1000. Some times some mate pairs are far off the actual value. This upper bound is used to get only good mate pairs. We may need to change the threshold for datasets with longer insert size.
									{
										insertSizes->push_back(locationOnEdgeRead1.at(k) - locationOnEdgeRead2.at(l));	// Insert the distance between the  reads in the list.
									}
							}
						}
					}
				}
			}
		}

		INT64 sum = 0, variance=0;
		if(insertSizes->size() == 0) // If no insert size found
		{
			cout << "No insert-size found for dataset: " << d << endl;
			meanOfInsertSizes.push_back(0);
			sdOfInsertSizes.push_back(0);
			continue;
		}

		for(UINT64 i = 0; i < insertSizes->size(); i++)
			sum += insertSizes->at(i);			// Calculate the sum of all the insert sizes of the current dataset.

		meanOfInsertSizes.push_back(sum/insertSizes->size());	// Calculate the mean. In push it in the variable of the OverlapGraph object.

		for(UINT64 i = 0; i < insertSizes->size(); i++)	// Calculate the variance of the current dataset.
				variance += (meanOfInsertSizes.at(d) - insertSizes->at(i)) * (meanOfInsertSizes.at(d) - insertSizes->at(i));

		sdOfInsertSizes.push_back(sqrt(variance/insertSizes->size()));  // Calculate and insert the standard deviation.

		// Print the values of the current dataset.
		cout << "Mean set to: " << meanOfInsertSizes.at(d) << endl;
		cout << "SD set to: " << sdOfInsertSizes.at(d) << endl;
		cout << "Reads on same edge: " << insertSizes->size() << endl;

	}
	for(UINT64 i = 0; i < meanOfInsertSizes.size(); i++)
	{
		if(longestMeanOfInsertSize < meanOfInsertSizes.at(i))
		{
			longestMeanOfInsertSize = meanOfInsertSizes.at(i);
		}
	}
	cout << "Mean of longest insert size: " << longestMeanOfInsertSize << endl;

	delete insertSizes;
	CLOCKSTOP;
	return true;
}

// CP: saveGraphToFile and readGraphFromFile are a pair of functions used together

/**********************************************************************************************************************
	Save the unitig graph in a text file
	Only stores the edges, there are two copies of each edge (forward and reverse). Only store one copy from smaller ID to larger ID here. when retrieve information, make two copies
	Sample:
12				// source read
14115			// destination read
2				// orientation of the edge
84				// overlap offset from the beginning of the source read to the beginning of the destination read
11				// number of read in the edge
27591			// the first read in the edge
7				// overlap offset from the source read to the first read (beginning to beginning)
0				// the orientation of the read
27049			// the second read in the edge
1				// overlap offset from the first read to the second read
1				// the orientation
// this will continue for all 11 reads in this edge and then start the next edge
**********************************************************************************************************************/
/*bool OverlapGraph::saveGraphToFile(string fileName)
{
	CLOCKSTART;
	ofstream filePointer;
	filePointer.open(fileName.c_str());
	if(filePointer == NULL)
		MYEXIT("Unable to open file: "+fileName);

	vector<UINT64> *list = new vector<UINT64>;
	for(UINT64 i = 1; i < graph->size(); i++) // for each node
	{
		if(!graph->at(i)->empty())
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)	// for each edge of the node
			{
				Edge * e = graph->at(i)->at(j);
				UINT64 source = e->getSourceRead()->getReadNumber();
				UINT64 destination = e->getDestinationRead()->getReadNumber();
				if(source < destination || (source == destination && e < e->getReverseEdge()))
				{
					list->push_back(source);	// store the edge information first
					list->push_back(destination);
					list->push_back(e->getOrientation());
					list->push_back(e->getOverlapOffset());
					list->push_back(e->getListOfReads()->size());	// store the size of the vector of elements in the edge.
					for(UINT64 k = 0; k < e->getListOfReads()->size(); k++)	// store information about the reads in the edge.
					{
						list->push_back(e->getListOfReads()->at(k));
						list->push_back(e->getListOfOverlapOffsets()->at(k));
						list->push_back(e->getListOfOrientations()->at(k));
					}
				}
			}
		}
	}

	for(UINT64 i = 0; i < list->size(); i++)	// store in a file for future use.
		filePointer<<list->at(i)<<endl;
	filePointer.close();
	delete list;
	CLOCKSTOP;
	return true;
}*/


bool OverlapGraph::saveGraphToFile(string fileName)
{
	CLOCKSTART;
	ofstream filePointer;
	filePointer.open(fileName.c_str());
	if(filePointer == NULL)
		MYEXIT("Unable to open file: "+fileName);
	filePointer << dataSet->getNumberOfUniqueReads() << endl;	// save the number of unique reads.
	filePointer << (int)(flowComputed) << endl;					// save whether the flow is computed or not
	filePointer << dataSet->numberOfPairedDatasets<< endl;		// number of pair end datasets.
	filePointer << meanOfInsertSizes.size()<< endl;				// Number of insert sizes for different datasets
	for(UINT64 i = 0; i< meanOfInsertSizes.size();i++)			// save the mean and sd of insert sizes if already calculated.
	{
		filePointer << meanOfInsertSizes.at(i) << endl;
		filePointer << sdOfInsertSizes.at(i) << endl;

	}
	for(UINT64 i = 1; i < graph->size(); i++)					// store the reads.
	{
		Read *r = dataSet->getReadFromID(i);					// get each read from the datasset
		filePointer << r->getReadNumber() << endl;				// save the read ID
		filePointer << r->getStringForward() << endl;			// save the forward string in the read.
		filePointer << r->superReadID << endl;					// save the super read ID.
		filePointer << r->getFrequency() << endl;				// save the frequecny of the read.
		filePointer << r->getMatePairList()->size() << endl;	// How many matepairs of the current read.
		for(UINT64 j=0; j< r->getMatePairList()->size(); j++ )	// store each of the matepair informatoin.
		{
			filePointer << r->getMatePairList()->at(j).matePairID << endl;										// ID of the matepair
			filePointer << (int)(r->getMatePairList()->at(j).matePairOrientation) << endl;						// Orientation of the matepair
			filePointer << (int)(r->getMatePairList()->at(j).datasetNumber) << endl;							// dataset number of the matepair
		}
	}
	vector<UINT64> *list = new vector<UINT64>;		// Now store the edges.
	for(UINT64 i = 1; i < graph->size(); i++) // for each node
	{
		if(!graph->at(i)->empty())
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)	// for each edge of the node
			{
				Edge * e = graph->at(i)->at(j);
				UINT64 source = e->getSourceRead()->getReadNumber();			// source read.
				UINT64 destination = e->getDestinationRead()->getReadNumber();	// destination read.
				if(source < destination || (source == destination && e < e->getReverseEdge()))
				{
					list->push_back(source);	// store the source read ID
					list->push_back(destination);	// destination read ID
					list->push_back(e->getOrientation());	// Orientaiton of the edge.
					list->push_back(e->getOverlapOffset());		// overlap offset of the edge
					list->push_back(e->flow);					// flow of the edge
					list->push_back(e->getListOfReads()->size());	// store the size of the vector of elements in the edge.
					for(UINT64 k = 0; k < e->getListOfReads()->size(); k++)	// store information about the reads in the edge.
					{
						list->push_back(e->getListOfReads()->at(k));				// store the ordered reads in the edge.
						list->push_back(e->getListOfOverlapOffsets()->at(k));		// ofset of this read from the previous read.
						list->push_back(e->getListOfOrientations()->at(k));			// orientiaon of the read in the edge.
					}
				}
			}
		}
	}

	for(UINT64 i = 0; i < list->size(); i++)	// store in a file for future use.
		filePointer<<list->at(i)<<endl;
	filePointer.close();
	delete list;
	CLOCKSTOP;
	return true;
}





/**********************************************************************************************************************
	Read the unitig graph from a binary file
**********************************************************************************************************************/
/*bool OverlapGraph::readGraphFromFile(string fileName)
{
	CLOCKSTART;
	ifstream filePointer;
	filePointer.open(fileName.c_str());
	if(filePointer == NULL)
		MYEXIT("Unable to open file: "+fileName);

	graph = new vector< vector<Edge *> * >;
	// Ted: Because read ID starts with 1, not 0, reserve numberOfUniqueReads + 1
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);

	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // initialize the list
	{
		vector<Edge *> *newList = new vector<Edge *>;
		graph->push_back(newList);
	}

	vector<UINT64> *list =  new vector<UINT64>;

	while(filePointer.good())	// read from the file.
	{
		UINT64 temp;
		filePointer>>temp;
		list->push_back(temp);
	};
	filePointer.close();

	for(UINT64 i = 0; i < list->size() -1;) // parse the numbers.
	// Ted: above loop iterate from 0 to size()-2, so size()-1 times..Why?
	{

		// Ted: The unitig grap should be 5 columns text file
		// Column 1: source read id.
		// Column 2:
		UINT64 source = list->at(i++);	// first number is the source read id.
		UINT64 destination = list->at(i++);	// destination read id.
		UINT64 orientation = list->at(i++);	// orientation of the edge.
		UINT64 overlapOffset = list->at(i++);	// overlap offset of the edge
		UINT64 numberOfReadsInEdge = list->at(i++);	// Array size
		vector<UINT64> *listReads = new vector<UINT64>;
		vector<UINT16> *listOverlapOffsets = new vector<UINT16>;
		vector<UINT8> *listOrientations = new vector<UINT8>;
		UINT64 length = 0;
		for(UINT64 j = 0; j < 3 * numberOfReadsInEdge; j += 3)	// Read all the three arrays. list of reads, overlap offsets and orientation of the reads.
		{
			listReads->push_back(list->at(i + j));
			listOverlapOffsets->push_back(list->at(i + j + 1));
			listOrientations->push_back(list->at(i + j + 2));
			length += list->at(i + j + 1);
		}

		Edge *edgeForward = new Edge();	// create an edge
		Edge *edgeReverse = new Edge();
		Read *read1, *read2;

		read1 = dataSet->getReadFromID(source);
		read2 = dataSet->getReadFromID(destination);


		vector<UINT64> *listReadsReverse = new vector<UINT64>;
		vector<UINT16> *listOverlapOffsetsReverse = new vector<UINT16>;
		vector<UINT8> *listOrientationsReverse = new vector<UINT8>;
		UINT64 size = listReads->size();

		UINT64 length1, length2, overlapOffsetForward, revereseOverlap;

		for(UINT64 j = 0; j < size; j++)	// creat the list of reads, overlap offsets and orientations for the reverse edge based on the number in the forward edge.
		{
			listReadsReverse->push_back(listReads->at(size-j-1));

			if(j == 0) // last/first read
			{
				length1 = read2->getReadLength();
				overlapOffsetForward = overlapOffset - length;
			}
			else // any read within the edge
			{
				length1 = dataSet->getReadFromID(listReads->at(size-j))->getReadLength();
				overlapOffsetForward = listOverlapOffsets->at(size-j);
			}
			length2 = dataSet->getReadFromID(listReads->at(size-j-1))->getReadLength();
			revereseOverlap = length1 + overlapOffsetForward - length2;
			listOverlapOffsetsReverse->push_back(revereseOverlap);

			listOrientationsReverse->push_back(!(listOrientations->at(size-j-1)));
		}
		UINT64 revereseOverlapOffset =  overlapOffset + read2->getReadLength() - read1->getReadLength();

		edgeForward->makeEdge(read1, read2, orientation, overlapOffset, listReads, listOverlapOffsets, listOrientations); // make the forward edge.
		edgeReverse->makeEdge(read2, read1, twinEdgeOrientation(orientation), revereseOverlapOffset, listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);	// make the reverse edge.

		edgeForward->setReverseEdge(edgeReverse); // set the reverse edge pointer
		edgeReverse->setReverseEdge(edgeForward); // set the reverse edge pinter.

		insertEdge(edgeForward);	// insert the forward edge
		insertEdge(edgeReverse);	// insert the reverse edge.


		i += numberOfReadsInEdge * 3;
	}
	delete list;
	CLOCKSTOP;
	return true;
}*/

bool OverlapGraph::readGraphFromFile(string fileName)
{
	CLOCKSTART;
	ifstream filePointer;
	filePointer.open(fileName.c_str());
	if(filePointer == NULL)
		MYEXIT("Unable to open file: "+fileName);


	UINT64 numOfUniqueRead,numOfMatePairDataset,temp1,temp2;
	filePointer >> numOfUniqueRead;							// Number of unique reads in the dataset.
	dataSet->numberOfUniqueReads = numOfUniqueRead;
	filePointer >> temp1;									// Flag for flow computed or not.
	flowComputed = temp1;
	filePointer >> numOfMatePairDataset;					// Number of mate pair datasets
	dataSet->numberOfPairedDatasets = numOfMatePairDataset;
	filePointer>>temp1;										// Number of mean and sd of insert size calculated.
	for(UINT64 i =0; i<temp1;i++)
	{
		filePointer>> temp2;								// get the mean of each dataset if calculated already.
		this->meanOfInsertSizes.push_back(temp2);
		filePointer>> temp2;								// get the standard deviation of insert size.
		this->sdOfInsertSizes.push_back(temp2);
	}
	for(UINT64 i = 0; i < numOfUniqueRead; i++)				// we have to read all the unique reads from the file.
	{
		UINT64 n;
		string s;
		Read *r = new Read;

		filePointer >> n;									// first number is the read ID
		r->setReadNumber(n);
		filePointer >> s;									// second string is the forward string of the read.
		r->setRead(s);
		filePointer >> n;									// thirs number is the ID of superRead. 0 means no super read.
		r->superReadID = n;
		filePointer >> n;									// frequency of this read in all the datasets.
		r->setFrequency(n);
		filePointer >> n;									// Number of mate pairs of this read.
		for(UINT64 j=0; j< n; j++ )							// Read all the matepair information.
		{
			UINT64 matePairID,orient,datasetNumber;
			filePointer >> matePairID;						// This is the ID of the matepair.
			filePointer >> orient;							// Orientation of the matepair.
			filePointer >> datasetNumber;					// dataset number of the matepair.
			r->addMatePair(matePairID,orient,datasetNumber);		// add the mate pair
		}
		dataSet->addRead(r);								// add the read in the dataset
	}
	graph = new vector< vector<Edge *> * >;
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);
	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // initialize the list
	{
		vector<Edge *> *newList = new vector<Edge *>;
		graph->push_back(newList);
	}
	vector<UINT64> *list =  new vector<UINT64>;
	while(filePointer.good())	// read from the file.
	{
		UINT64 temp;
		filePointer>>temp;
		list->push_back(temp);
	};
	filePointer.close();
	for(UINT64 i = 0; i < list->size() -1;) // parse the numbers.
	{

		UINT64 source = list->at(i++);	// first number is the source read id.
		UINT64 destination = list->at(i++);	// destination read id.
		UINT64 orientation = list->at(i++);	// orientation of the edge.
		UINT64 overlapOffset = list->at(i++);	// overlap offset of the edge
		UINT64 flow = list->at(i++); // get the flow
		UINT64 numberOfReadsInEdge = list->at(i++);	// Array size
		vector<UINT64> *listReads = new vector<UINT64>;
		vector<UINT16> *listOverlapOffsets = new vector<UINT16>;
		vector<UINT8> *listOrientations = new vector<UINT8>;
		UINT64 length = 0;
		for(UINT64 j = 0; j < 3 * numberOfReadsInEdge; j += 3)	// Read all the three arrays. list of reads, overlap offsets and orientation of the reads.
		{
			listReads->push_back(list->at(i + j));					// This is the read ID
			listOverlapOffsets->push_back(list->at(i + j + 1));		// This is the offset of the read from the previous read.
			listOrientations->push_back(list->at(i + j + 2));		// This is the orientaion of the read on this edge.
			length += list->at(i + j + 1);
		}

		Edge *edgeForward = new Edge();	// create an edge
		Edge *edgeReverse = new Edge();
		Read *read1, *read2;

		read1 = dataSet->getReadFromID(source);
		read2 = dataSet->getReadFromID(destination);


		vector<UINT64> *listReadsReverse = new vector<UINT64>;
		vector<UINT16> *listOverlapOffsetsReverse = new vector<UINT16>;
		vector<UINT8> *listOrientationsReverse = new vector<UINT8>;
		UINT64 size = listReads->size();

		UINT64 length1, length2, overlapOffsetForward, revereseOverlap;

		for(UINT64 j = 0; j < size; j++)	// creat the list of reads, overlap offsets and orientations for the reverse edge based on the number in the forward edge.
		{
			listReadsReverse->push_back(listReads->at(size-j-1));

			if(j == 0) // last/first read
			{
				length1 = read2->getReadLength();
				overlapOffsetForward = overlapOffset - length;
			}
			else // any read within the edge
			{
				length1 = dataSet->getReadFromID(listReads->at(size-j))->getReadLength();
				overlapOffsetForward = listOverlapOffsets->at(size-j);
			}
			length2 = dataSet->getReadFromID(listReads->at(size-j-1))->getReadLength();
			revereseOverlap = length1 + overlapOffsetForward - length2;
			listOverlapOffsetsReverse->push_back(revereseOverlap);

			listOrientationsReverse->push_back(!(listOrientations->at(size-j-1)));
		}
		UINT64 revereseOverlapOffset =  overlapOffset + read2->getReadLength() - read1->getReadLength();

		edgeForward->makeEdge(read1, read2, orientation, overlapOffset, listReads, listOverlapOffsets, listOrientations); // make the forward edge.
		edgeReverse->makeEdge(read2, read1, twinEdgeOrientation(orientation), revereseOverlapOffset, listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);	// make the reverse edge.

		edgeForward->setReverseEdge(edgeReverse); // set the reverse edge pointer
		edgeForward->flow = flow;				// Add the flow in the forward edge.
		edgeReverse->setReverseEdge(edgeForward); // set the reverse edge pinter.
		edgeReverse->flow = flow;					// add the flow in the reverse edge.

		insertEdge(edgeForward);	// insert the forward edge
		insertEdge(edgeReverse);	// insert the reverse edge.


		i += numberOfReadsInEdge * 3;
	}
	delete list;
	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
	Calculate min cost flow

	// An sample input to the CS2 algorithm

	p min       3840      13449										// p min numberOfNodes numberOfEdges
	n          1         0											// initial flow of 0 in node 1, node 1 is the supersource
	n       3840         0											// initial flow of 0 in node 3840, which is the supersink (equal to the number of nodes)
	a       3840          1          1    1000000    1000000		// edge from supersink to supersource, LowBound(1), UpperBound(1000000), Cost per unit of flow(1000000)
	a          1          2          0    1000000          0		// edge from supersource to node 2, with the defined cost function
	// this continues to
	// connect each node to supersource and supersink
	// connect every edge in the original graph

**********************************************************************************************************************/
bool OverlapGraph::calculateFlow(string inputFileName, string outputFileName)
{
	CLOCKSTART;
	for(UINT64 i = 1; i < graph->size(); i++) // clear all the flow
	{
		if(!graph->at(i)->empty())
		{
			for(UINT64 j = 0; j< graph->at(i)->size(); j++)
			{
				graph->at(i)->at(j)->flow =0;
			}

		}
	}
	// CP: the number of vertices for CS2 graph, two CS2 nodes for each read nodes plus supersource and supersink
	UINT64 V = numberOfNodes * 2 + 2;
	// CP: the number of edges for CS2 graph,
	// CP: 3 CS2 edges for each overlap graph edge, plus two edges for CS nodes to supersource and supersink, plus an edge between supersource and supersink
	UINT64 E = numberOfEdges * 3 + numberOfNodes * 4 + 1;
	UINT64 SUPERSOURCE = 1;
	UINT64 SUPERSINK = V;
	// For each edge in the overlap graph, we add 3 edges in the directed graph. Two nodes are created for each node in the original graph.
	// A super source and a super sink is added. Each node is connected to the super source and super sink.
	INT64 FLOWLB[3], FLOWUB[3], COST[3];			// Flow bounds and cost of the edges.
	// CP: comment on what's inputfile and output files.
	// CP: it's confusing to open outputFile with inputFileName
	// BH: This is an ouput file for this function, but an input file for CS2
	ofstream outputFile;
	outputFile.open(inputFileName.c_str());
	if(outputFile == NULL)
		MYEXIT("Unable to open file: "+inputFileName);

	stringstream ss;
	ss << "p min " << setw(10) << V << " " << setw(10) << E << endl;  	// Number of nodes and edges in the new graph.
	ss << "n " << setw(10) << SUPERSOURCE << setw(10) << " 0" << endl;	// Flow in the super source
	ss << "n " << setw(10) << SUPERSINK << setw(10) << " 0" << endl;	// Flow in the super sink.
	FLOWLB[0] = 1; FLOWUB[0] = 1000000; COST[0] = 1000000;
	ss << "a " << setw(10) << SUPERSINK << " " << setw(10) << SUPERSOURCE << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl; // Add an edge from super sink to super source with very high cost.

	// CP: this is a lookup table from the overlap graph node ID to the CS2 graph ID
	vector<UINT64> *listOfNodes = new vector<UINT64>;
	// CP: this is a lookup table from the CS2 graph ID to the overlap graph node ID
	vector<UINT64> *listOfNodesReverse = new vector<UINT64>;

	// For n nodes in the graph, CS2 requires that the nodes are numbered from 1 to n. In the overlap graph, the nodes does not have sequencinal ID. We need to convert them to 1 - n

	// If the ID of a node in the original graph is 100 and directed graph is 5
	// Then listOfNodes->at(100) is equal to 5
	// and ListOfNodesReverse->at(5) is equal to 100.

	for(UINT64 i = 0; i <= graph->size(); i++)
	{
		listOfNodes->push_back(0);
		listOfNodesReverse->push_back(0);
	}

	// This loop set lower bound and upper bound of each node to super source and to super sink. All costs are 0.
	UINT64 currentIndex = 1;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			listOfNodes->at(i) = currentIndex;					// Mapping between original node ID and cs2 node ID
			listOfNodesReverse->at(currentIndex) = i;			// Mapping between original node ID and cs2 node ID
			// CP: don't you create two CS2 node for each original node? where is the second CS2 node ID?
			// BH: Yes I created two nodes for CS2. For a node u in the listOfNodes. We created 2*u and 2*u+1 in for the cs2
			FLOWLB[0] = 0; FLOWUB[0] = 1000000; COST[0] = 0;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << 2 * currentIndex << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			currentIndex++;
		}
	}

	// This loop converts the original bi-directed edges to directed edges (1 becomes 6).
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{

				// CP: u and v are the CS2 node IDs of the source node and destination node, respectively, of the edge
				Edge *edge = graph->at(i)->at(j);
				UINT64 u = listOfNodes->at(edge->getSourceRead()->getReadNumber());
				UINT64 v = listOfNodes->at(edge->getDestinationRead()->getReadNumber());

				// set the bound and cost here
				// if edge has more than 20 reads:
				//   FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
				//   FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
				//   FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
				// else:
				//   FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
				//   FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
				//   FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
				// cost function is set in such a way that for the composite edges with more that 20 reads in them we will push at least one units of flow.
				// Otherwise we set the lower bound of flow to 0, meaning that these edges might not have any flow at the endl.
				// The cost of pushing the first using of flow in very cheap and then we pay high price to push more than 1 units of flow.
				// This will ensure that we do not push flow were it is not necessary.
				// CP: if we need to change the cost function, we just need to change this function, right?
				// BH: Yes, we only need to change this function if we want to use different cost function.
				calculateBoundAndCost(edge, FLOWLB, FLOWUB, COST);


				if(u < v || (u == v && edge < edge->getReverseEdge()))
				{
					// Here for each edge we add three edges with different values of cost and bounds.
					// Total 6 edges considering the reverse edge too.
					// For details on how to convert the edges off different types please see my thesis.

					// BH: u1, u2, v1 and v2 are the actual CS2 node IDs
					UINT64 u1 = 2 * u, u2 = 2 * u + 1, v1 =  2 * v, v2 = 2 * v + 1;

					// Connect the edges of the original graph
					// for each orignal edge, we add six edges
					if(edge->getOrientation() == 0)
					{
						// first edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						// second edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						// third edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
					}
					else if(edge->getOrientation() == 1)
					{
						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
					else if(edge->getOrientation() == 2)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
					else if(edge->getOrientation() == 3)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
				}
			}
		}
	}
	outputFile << ss.str();		// Write the string in a file for CS2
	outputFile.close();

	// CP: clear ss, right?
	ss.str(std::string());

	// CP: what are  you doing here?
	// BH: CS2 requires the file name to be char * not string. We convert the string filename to char *
	char * inFile = new char[inputFileName.size() + 1];
	std::copy(inputFileName.begin(), inputFileName.end(), inFile);
	inFile[inputFileName.size()] = '\0';

	char * outFile = new char[outputFileName.size() + 1];
	std::copy(outputFileName.begin(), outputFileName.end(), outFile);
	outFile[outputFileName.size()] = '\0';


	cout << "Calling CS2" << endl;
	main_cs2(inFile,outFile);			// Call CS2
	cout << "CS2 finished" << endl;

	delete[] inFile;
	delete[] outFile;

	ifstream inputFile;
	inputFile.open(outputFileName.c_str());
	if(inputFile == NULL)
		MYEXIT("Unable to open file: "+outputFileName);

	UINT64 lineNum = 0;
	while(!inputFile.eof())
	{
		lineNum ++;
		UINT64 source, destination, flow;
		inputFile >> source >> destination >> flow;		// get the flow from CS2
		// CP: give an sample of the CS2 output
		/*
		From To Flow
		1 2421 0 from node 1 to node 2421 with flow of 0
		1 3 0	from node 1 to node 3 with flow of 0
		*/
		if(source != SUPERSINK && source != SUPERSOURCE && destination != SUPERSOURCE && destination != SUPERSINK && flow!=0)
		{
			UINT64 mySource = listOfNodesReverse->at(source/2);				// Map the source to the original graph
			UINT64 myDestination = listOfNodesReverse->at(destination/2);	// Map the destination in the original graph
			Edge *edge = findEdge(mySource, myDestination);					// Find the edge in the original graph.
			edge->flow += flow;												// Add the flow in the original graph.
		}
	}
	inputFile.close();
	delete listOfNodes;
	delete listOfNodesReverse;
	this->flowComputed = true;
	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
	This function calculates the cost and bounds for an edge in the overlap graph.
	This function is very sensitive to the assembled contigs and to the cost function parameters
	BH: changing the bounds and threshold of number of nodes in the edge (here 20) many give use very wrong flow.

	CP: given an *edge, calculate and return FLOWLB, FLOWUB, and COST
	CP: FLOWLB, FLOWUB, and COST are all array of size 3 because each overlap graph edge is represented by 3 CS2 edges to define a cost function
**********************************************************************************************************************/
bool OverlapGraph::calculateBoundAndCost(const Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST)
{
	for(UINT64 i = 0; i < 3; i++)		// For the simple edges we put very high cost
	{
		FLOWLB[i] = 0; FLOWUB[i] = 10; COST[i] = 500000;
	}

	if(!edge->getListOfReads()->empty()) // Composite Edge
	{
		if(edge->getListOfReads()->size() > 20) // Composite Edge of at least 20 reads. Must have at least one unit of flow.
		{
			// the first 1 flow must be used and has a cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carry the first 1 flow
			FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second 1 flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the second flow
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
		else // Short composite edge containing less than 20 reads. May have zero flow.
		{
			// the first 1 flow may not be required, but has a low cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carries the first 1 flow
			FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second unit of flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the two units of flow.
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
	}

	return true;
}


/**********************************************************************************************************************
	Calculate min cost flow

	// An sample input to the CS2 algorithm

	p min       3840      13449										// p min numberOfNodes numberOfEdges
	n          1         0											// initial flow of 0 in node 1, node 1 is the supersource
	n       3840         0											// initial flow of 0 in node 3840, which is the supersink (equal to the number of nodes)
	a       3840          1          1    1000000    1000000		// edge from supersink to supersource, LowBound(1), UpperBound(1000000), Cost per unit of flow(1000000)
	a          1          2          0    1000000          0		// edge from supersource to node 2, with the defined cost function
	// this continues to
	// connect each node to supersource and supersink
	// connect every edge in the original graph

**********************************************************************************************************************/
bool OverlapGraph::calculateFlow2(string inputFileName, string outputFileName)
{
	CLOCKSTART;
	for(UINT64 i = 1; i < graph->size(); i++) // clear all the flow
	{
		if(!graph->at(i)->empty())
		{
			for(UINT64 j = 0; j< graph->at(i)->size(); j++)
			{
				graph->at(i)->at(j)->flow =0;
			}

		}
	}
	// CP: the number of vertices for CS2 graph, two CS2 nodes for each read nodes plus supersource and supersink
	UINT64 V = numberOfNodes * 2 + 2;
	// CP: the number of edges for CS2 graph,
	// CP: 3 CS2 edges for each overlap graph edge, plus two edges for CS nodes to supersource and supersink, plus an edge between supersource and supersink
	UINT64 E = numberOfEdges * 3 + numberOfNodes * 4 + 1;
	UINT64 SUPERSOURCE = 1;
	UINT64 SUPERSINK = V;
	// For each edge in the overlap graph, we add 3 edges in the directed graph. Two nodes are created for each node in the original graph.
	// A super source and a super sink is added. Each node is connected to the super source and super sink.
	INT64 FLOWLB[3], FLOWUB[3], COST[3];			// Flow bounds and cost of the edges.
	// CP: comment on what's inputfile and output files.
	// CP: it's confusing to open outputFile with inputFileName
	// BH: This is an ouput file for this function, but an input file for CS2
	ofstream outputFile;
	outputFile.open(inputFileName.c_str());
	if(outputFile == NULL)
		MYEXIT("Unable to open file: "+inputFileName);

	stringstream ss;
	ss << "p min " << setw(10) << V << " " << setw(10) << E << endl;  	// Number of nodes and edges in the new graph.
	ss << "n " << setw(10) << SUPERSOURCE << setw(10) << " 0" << endl;	// Flow in the super source
	ss << "n " << setw(10) << SUPERSINK << setw(10) << " 0" << endl;	// Flow in the super sink.
	FLOWLB[0] = 1; FLOWUB[0] = 1000000; COST[0] = 1000000;				// Cost of unit of flow from super sink to super source
	ss << "a " << setw(10) << SUPERSINK << " " << setw(10) << SUPERSOURCE << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl; // Add an edge from super sink to super source with very high cost.

	// CP: this is a lookup table from the overlap graph node ID to the CS2 graph ID
	vector<UINT64> *listOfNodes = new vector<UINT64>;
	// CP: this is a lookup table from the CS2 graph ID to the overlap graph node ID
	vector<UINT64> *listOfNodesReverse = new vector<UINT64>;

	// For n nodes in the graph, CS2 requires that the nodes are numbered from 1 to n. In the overlap graph, the nodes does not have sequencinal ID. We need to convert them to 1 - n

	// If the ID of a node in the original graph is 100 and directed graph is 5
	// Then listOfNodes->at(100) is equal to 5
	// and ListOfNodesReverse->at(5) is equal to 100.

	for(UINT64 i = 0; i <= graph->size(); i++)
	{
		listOfNodes->push_back(0);
		listOfNodesReverse->push_back(0);
	}
	//markEdgeThatWillHaveOneFlow();
	// This loop set lower bound and upper bound of each node to super source and to super sink. All costs are 0.
	UINT64 currentIndex = 1;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			listOfNodes->at(i) = currentIndex;					// Mapping between original node ID and cs2 node ID
			listOfNodesReverse->at(currentIndex) = i;			// Mapping between original node ID and cs2 node ID
			// CP: don't you create two CS2 node for each original node? where is the second CS2 node ID?
			// BH: Yes I created two nodes for CS2. For a node u in the listOfNodes. We created 2*u and 2*u+1 in for the cs2
			FLOWLB[0] = 0; FLOWUB[0] = 1000000; COST[0] = 0;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << 2 * currentIndex << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			currentIndex++;
		}
	}

	// This loop converts the original bi-directed edges to directed edges (1 becomes 6).
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{

				// CP: u and v are the CS2 node IDs of the source node and destination node, respectively, of the edge
				Edge *edge = graph->at(i)->at(j);
				UINT64 u = listOfNodes->at(edge->getSourceRead()->getReadNumber());
				UINT64 v = listOfNodes->at(edge->getDestinationRead()->getReadNumber());

				// set the bound and cost here
				// if edge has more than 20 reads:
				//   FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
				//   FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
				//   FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
				// else:
				//   FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
				//   FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
				//   FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
				// cost function is set in such a way that for the composite edges with more that 20 reads in them we will push at least one units of flow.
				// Otherwise we set the lower bound of flow to 0, meaning that these edges might not have any flow at the endl.
				// The cost of pushing the first using of flow in very cheap and then we pay high price to push more than 1 units of flow.
				// This will ensure that we do not push flow were it is not necessary.
				// CP: if we need to change the cost function, we just need to change this function, right?
				// BH: Yes, we only need to change this function if we want to use different cost function.
				calculateBoundAndCost2(edge, FLOWLB, FLOWUB, COST);


				if(u < v || (u == v && edge < edge->getReverseEdge()))
				{
					// Here for each edge we add three edges with different values of cost and bounds.
					// Total 6 edges considering the reverse edge too.
					// For details on how to convert the edges off different types please see my thesis.

					// BH: u1, u2, v1 and v2 are the actual CS2 node IDs
					UINT64 u1 = 2 * u, u2 = 2 * u + 1, v1 =  2 * v, v2 = 2 * v + 1;

					// Connect the edges of the original graph
					// for each orignal edge, we add six edges
					if(edge->getOrientation() == 0)
					{
						// first edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						// second edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						// third edge in the cost function, forward and reverse
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
					}
					else if(edge->getOrientation() == 1)
					{
						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
					else if(edge->getOrientation() == 2)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
					else if(edge->getOrientation() == 3)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
				}
			}
		}
	}
	outputFile << ss.str();		// Write the string in a file for CS2
	outputFile.close();

	// CP: clear ss, right?
	ss.str(std::string());

	// CP: what are  you doing here?
	// BH: CS2 requires the file name to be char * not string. We convert the string filename to char *
	char * inFile = new char[inputFileName.size() + 1];
	std::copy(inputFileName.begin(), inputFileName.end(), inFile);
	inFile[inputFileName.size()] = '\0';

	char * outFile = new char[outputFileName.size() + 1];
	std::copy(outputFileName.begin(), outputFileName.end(), outFile);
	outFile[outputFileName.size()] = '\0';


	cout << "Calling CS2" << endl;
	main_cs2(inFile,outFile);			// Call CS2
	cout << "CS2 finished" << endl;

	delete[] inFile;
	delete[] outFile;

	ifstream inputFile;
	inputFile.open(outputFileName.c_str());
	if(inputFile == NULL)
		MYEXIT("Unable to open file: "+outputFileName);

	UINT64 lineNum = 0;
	while(!inputFile.eof())
	{
		lineNum ++;
		UINT64 source, destination, flow;
		inputFile >> source >> destination >> flow;		// get the flow from CS2
		// CP: give an sample of the CS2 output
		/*
		From To Flow
		1 2421 0 from node 1 to node 2421 with flow of 0
		1 3 0	from node 1 to node 3 with flow of 0
		*/
		if(source != SUPERSINK && source != SUPERSOURCE && destination != SUPERSOURCE && destination != SUPERSINK && flow!=0)
		{
			UINT64 mySource = listOfNodesReverse->at(source/2);				// Map the source to the original graph
			UINT64 myDestination = listOfNodesReverse->at(destination/2);	// Map the destination in the original graph
			Edge *edge = findEdge(mySource, myDestination);					// Find the edge in the original graph.
			edge->flow += flow;												// Add the flow in the original graph.
		}
	}
	inputFile.close();
	delete listOfNodes;
	delete listOfNodesReverse;
	this->flowComputed = true;
	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
	This function calculates the cost and bounds for an edge in the overlap graph.
	This function is very sensitive to the assembled contigs. CP: what does this mean?
	BH: changing the bounds and threshold of number of nodes in the edge (here 20) many give use very wrong flow.

	CP: given an *edge, calculate and return FLOWLB, FLOWUB, and COST
	CP: FLOWLB, FLOWUB, and COST are all array of size 3 because each overlap graph edge is represented by 3 CS2 edges to define a cost function
**********************************************************************************************************************/
bool OverlapGraph::calculateBoundAndCost2(Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST)
{
	for(UINT64 i = 0; i < 3; i++)		// For the simple edges we put very high cost
	{
		FLOWLB[i] = 0; FLOWUB[i] = 10; COST[i] = 500000;
	}

	getBaseByBaseCoverage(edge);
	if(!edge->getListOfReads()->empty()) // Composite Edge
	{
		//if(edge->hignCoverageAndMatepairFlag ==  true)
		if(edge->coverageDepth >=coverageDepthLB && edge->coverageDepth <=coverageDepthUB)
		{
			// the first 1 flow must be used and has a cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carry the first 1 flow
			FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
			// this edge carries the second 1 flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			// this edge provides additional flow after the second flow
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
		else if(edge->coverageDepth > coverageDepthUB)
		{
			// the first 1 flow may not be required, but has a low cost of 1
			// each additional flow up to 8 flows has a cost of 100000

			// this edge carries the first 1 flow
			FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1000 + ( edge->coverageDepth - coverageDepthUB ) * 1;
			// this edge carries the second unit of flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000 + ( edge->coverageDepth - coverageDepthUB ) * 50;
			// this edge provides additional flow after the two units of flow.
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000 + ( edge->coverageDepth - coverageDepthUB ) * 100;
		}
		else
		{
			FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1000 + ( coverageDepthLB - edge->coverageDepth  ) * 1;
			// this edge carries the second unit of flow with high cost
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000 + ( coverageDepthLB - edge->coverageDepth ) * 50;
			// this edge provides additional flow after the two units of flow.
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000 + ( coverageDepthLB - edge->coverageDepth ) * 100;

		}
	}

	return true;
}


/**********************************************************************************************************************
	This function will mark all the edges that has coverage depth within the interval [coverageDepthLB, coverageDepthUB]
	and also the edges that share matepair with these edges.
**********************************************************************************************************************/
void OverlapGraph::markEdgeThatWillHaveOneFlow()
{
	for(UINT64 i = 1; i < graph->size(); i++) // for each node in the graph
	{
		if(!graph->at(i)->empty()) // this node has some edges.
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++) // for each edge of the current node
			{
				Edge *edge1 = graph->at(i)->at(j);	// this is the current edge
				getBaseByBaseCoverage(edge1);		// calculate the coverage depth of the current edge
				if(edge1->coverageDepth >=coverageDepthLB && edge1->coverageDepth <=coverageDepthUB)	// if the coverage depth is within the range [coverageDepthLB, coverageDepthUB]
				{
					edge1->hignCoverageAndMatepairFlag = true; // Mark this edge. We will set lower bound of flow for these edges to 1
					for(UINT64 k = 0; k<edge1->getListOfReads()->size(); k++)	// Now talke all the reads in the edge
					{
						UINT64 readID = edge1->getListOfReads()->at(k);	// This is the reads ID
						Read *read = dataSet->getReadFromID(readID);	// get the read object from the dataset
						for(UINT64 l = 0; l < read->getMatePairList()->size();l++)	// Find all the matepairs of this read.
						{
							UINT64 matePairID = read->getMatePairList()->at(l).matePairID;	// Get the matepair ID
							Read *matePair = dataSet->getReadFromID(matePairID);			// Get the matepair read object.
							for(UINT64 m = 0; m < matePair->getListOfEdgesForward()->size(); m++)	// List of the edges that contains the matepair.
							{
									Edge *edge2 = matePair->getListOfEdgesForward()->at(m);	// Get the edge
									edge2->hignCoverageAndMatepairFlag = true;			// We also need to set the lower bound of flow to these edges to one since they share edge with edge that has our desired coverage depth
									edge2->getReverseEdge()->hignCoverageAndMatepairFlag = true; // Mark the reverse edge too. Will set the lower bound of flow to 1
							}
						}
					}
				}
			}
		}
	}
}





/**********************************************************************************************************************
	return edge between source and destination
**********************************************************************************************************************/
Edge * OverlapGraph::findEdge(UINT64 source, UINT64 destination)
{
	for(UINT64 i = 0; i < graph->at(source)->size(); i++) // For the list of edges of the source node.
	{
		if(graph->at(source)->at(i)->getDestinationRead()->getReadNumber() == destination)	// check if there is an edge to destination
			return graph->at(source)->at(i);	// return the edge.
	}
	cout << "Check for error " << source << " to " << destination << endl;
	MYEXIT("Unable to find edge");
}



/**********************************************************************************************************************
	Checks if there is an edge (source, destination)
**********************************************************************************************************************/
bool OverlapGraph::isEdgePresent(UINT64 source, UINT64 destination)
{
	for(UINT64 i = 0; i < graph->at(source)->size(); i++)	// flro the list of edges of the source node
	{
		if(graph->at(source)->at(i)->getDestinationRead()->getReadNumber() == destination)	// check if there is an edge to destination
			return true;	// return true if there is an edge (source,destination)
	}
	return false;	// edge not found between source and destination
}


/**********************************************************************************************************************
	This function finds all the paths between two reads in the overlap graph.
**********************************************************************************************************************/

// CP: inputs: read1 and read2 are the pair, orient of the pair is defined in the MPlist struct, datasetNumber retrievs the mean and SD of insert size of the dataset
// CP: outputs: copyOfPath is a path between the two reads and copyOfFlags indicates whether the connection between this edge to the next edge is supported by all possible paths or not
// CP: return true if valid paths are found between them.
bool OverlapGraph::findPathBetweenMatepairs(const Read * read1, const Read * read2, UINT8 orient, UINT8 datasetNumber, vector <Edge *> &copyOfPath, vector <UINT64> &copyOfFlags)
{
	UINT64 pathFound = 0;			// CP: the total number of paths found between the two reads
	// CP: two variables passed to the exploreGraph function for return values
	vector <Edge *> firstPath;
	vector <UINT64> flags;

	copyOfPath.clear();
	copyOfFlags.clear();
	//UINT64 flag = 0;

	vector<Edge *> *listRead1, *listRead2;
	vector<UINT64> *locationOnEdgeRead1, *locationOnEdgeRead2;

	// mate pair orientation
	// 0 = 00 means the reverse of this read and the reverse of the second read are matepairs.
	// 1 = 01 means the reverse of this read and the forward of the second read are matepairs.
	// 2 = 10 means the forward of this read and the reverse of the second read are matepairs.
	// 3 = 11 means the forward of this read and the forward of the second read are matepairs.

	listRead1 = (orient == 2 || orient == 3) ? read1->getListOfEdgesForward() : read1->getListOfEdgesReverse();
	locationOnEdgeRead1 = (orient == 2 || orient == 3) ? read1->getLocationOnEdgeForward() : read1->getLocationOnEdgeReverse();

	// CP: Explain how this forward-forward/forward-reverse is related to mate-pair orientation listed above?
	// BH: if the first read is a forward read and the second read is a reverse read then we have to look for the forward read in the first edge and the reverse read in the second edge.

	// CP: use this, if the matepairs are forward - forward
	// CP: Ted, we need to add this to the config file, instead of hard-coding it here
	//listRead2 = (orient == 1 || orient == 3) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse(); // if the matepairs are forward - forward
	//locationOnEdgeRead2 = (orient == 1 || orient == 3) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse(); // if the matepairs are forward - forward

	// CP: use this, if the matepairs are forward - reverse
	listRead2 = (orient == 0 || orient == 2) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse(); // if the matepairs are forward - reverse
	locationOnEdgeRead2 = (orient == 0 || orient == 2) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse(); // if the matepairs are forward - reverse

	// CP: return false, if the two reads are not part of any edge or they are on the same edge
	if(listRead1->size()==0 || listRead2->size()==0)	// Not  present in any edges.
	{
		return false;
	}
	else
	{
		for(UINT64 i = 0; i < listRead1->size(); i++)
		{
			for(UINT64 j = 0; j < listRead2->size(); j++)
			{
				Edge * firstEdge = listRead1->at(i);
				Edge * lastEdge = listRead2->at(j);
				if(firstEdge ==lastEdge || firstEdge == lastEdge->getReverseEdge()) // Both are on the same edge
				{
					return false;
				}
			}
		}
	}

	for(UINT64 i = 0; i < listRead1->size(); i++)	// If not on the same edge.
	{
		for(UINT64 j = 0; j < listRead2->size(); j++)
		{
			Edge * firstEdge = listRead1->at(i);
			Edge * lastEdge = listRead2->at(j);
			if(firstEdge !=lastEdge && firstEdge != lastEdge->getReverseEdge()) //not on the same edge
			{
				// distance from the end of the read1 to the end of the last read of the edge,
				UINT64 distanceOnFirstEdge = firstEdge->getOverlapOffset() - locationOnEdgeRead1->at(i);
				// distance from the beginning of the first read of the edge to the beginning of the read2
				UINT64 distanceOnLastEdge = locationOnEdgeRead2->at(j);
				// the two reads can't be too far apart within their own edges
				if(distanceOnFirstEdge + distanceOnLastEdge < getMean(datasetNumber) + insertSizeRangeSD * getSD(datasetNumber))
				{
					UINT64 newPaths= exploreGraph(firstEdge, lastEdge, distanceOnFirstEdge, distanceOnLastEdge, datasetNumber, 0, firstPath, flags);	// from firs edge  try to find a path to the last edge.
					if(newPaths > 0)	// If some path found.
					{
						pathFound+=newPaths;	// How many paths found.
						if(copyOfPath.empty())	// Path found for the first time.
						{
							for(UINT64 k = 0; k < firstPath.size(); k++)	// Add the paths in the list.
								copyOfPath.push_back(firstPath.at(k));
							for(UINT64 k = 0; k < firstPath.size() - 1; k++)	// Also add the flag.
								copyOfFlags.push_back(flags.at(k));
						}
						else		// Not the first path
						{
							// CP: compare each supported pair of edge in copyOfPath with all pairs in firstpath
							// CP: if this supported pair of edge in copyOfPath is still supported by a pair in firstPath, it remains supported.
							// CP: if this supported pair of edge in copyOfPath is not supported by any pair in firstPath, it is changed to not supported.
							UINT64 k , l;
							for( k = 0; k< copyOfPath.size() - 1; k++)	// Check if the previous path contains the same pair of edges adjacent to the new path and has flag 1.
							{
								for( l = 0; l < firstPath.size() - 1; l++)
								{
									if(copyOfPath.at(k) == firstPath.at(l) &&  copyOfPath.at(k+1) == firstPath.at(l+1) && flags.at(l) == 1)
										break;
								}
								if(l == firstPath.size() - 1)	// This pair is not supportd
									copyOfFlags.at(k) = 0;
							}
						}
					}
				}
			}
		}
	}
	return true;
}

/**********************************************************************************************************************
	This function returns the edit distance between two strings.
	Code downloaded from http://rosettacode.org/wiki/Levenshtein_distance#C.2B.2B
**********************************************************************************************************************/
UINT64 OverlapGraph::calculateEditDistance(const  string & s1, const string & s2)
{
	const UINT64 m(s1.size());
	const UINT64 n(s2.size());
	if( m==0 )
		return n;
	if( n==0 )
		return m;
	UINT64 *costs = new UINT64[n + 1];
	for( UINT64 k=0; k<=n; k++ )
		costs[k] = k;

	UINT64 i = 0;
	for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
	{
		costs[0] = i+1;
		UINT64 corner = i;
		UINT64 j = 0;
		for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
		{
			UINT64 upper = costs[j+1];
			if( *it1 == *it2 )
			{
				costs[j+1] = corner;
			}
			else
			{
				UINT64 t(upper<corner?upper:corner);
				costs[j+1] = (costs[j]<t?costs[j]:t)+1;
			}
			corner = upper;
		}
	}
	UINT64 result = costs[n];
	delete [] costs;
	//cout << s1 << endl << s2 << endl << result<< endl;
	return result;
}


/**********************************************************************************************************************
	Explore all paths starting from firstEdge and tries to reach lastEdge using depth first search.
	Depth is limited to 100 to reduce running time
**********************************************************************************************************************/
// CP: this a recursive function
// CP: inputs: firstEdge and lastEdge in the path, distanceOnFirstEdge and distanceOnLastEdge are the relative position of the paired reads on the two edges
// CP: inputs: datasetNumber retrieves the mean and SD of insert size of the dataset
// CP: outputs: level is the level of the depth first search, firstPath and flags are ???
// BH: when we find the first path we save it and flag is used to indicate that all the edges in the first path are supported.
//     When we find another path, we unmark the flag for the pair of edges in the first path that are not present in the second path and so on.
// CP: return the number of paths found
UINT64 OverlapGraph::exploreGraph(Edge* firstEdge, Edge * lastEdge, UINT64 distanceOnFirstEdge, UINT64 distanceOnLastEdge, UINT64 datasetNumber, UINT64 level, vector <Edge *> &firstPath, vector <UINT64> &flags)
{
	// CP: use static variables to carry their values through the recursive calls of this function
	static UINT64 pathFound;					// number of paths found
	static vector <Edge *> listOfEdges;			// list of edges in the current path
	static vector <UINT64> pathLengths;			// length of the path from the beginning of the source read of the first edge to the current edge's source read's beginning
	if(level == 0)
	{
		// clear the variables and resize their capacity to free memory
		pathFound = 0;
		firstPath.resize(0);
		flags.resize(0);
		listOfEdges.resize(0);
		pathLengths.resize(0);
	}
	else
	{
		listOfEdges.resize(level);
		pathLengths.resize(level);
	}

	// BH: we do not go deeper than 100 levels. We can put this in the config file.
	// CP: when reaching the maximum depth, return 0 path found and exit the recursive call loop
	if(level > 100) return 0; // Do not go very deep.


	if(level == 0)
	{
		// CP: at level 0, start from the end of the first edge
		listOfEdges.push_back(firstEdge);
		pathLengths.push_back(distanceOnFirstEdge);
	}
	else
	{
		if(firstEdge == lastEdge) // Destination found.
		{
			// If we read our destination read, we check if the distance is within 3 sd of the mean.
			// mean - 3*sd can be negative. I think we do not have to worry about it.
			if((INT64)(distanceOnLastEdge + pathLengths.at(level - 1)) >= (INT64)((INT64)(getMean(datasetNumber)) - insertSizeRangeSD * (INT64)(getSD(datasetNumber))) && distanceOnLastEdge + pathLengths.at(level - 1) <= getMean(datasetNumber) + insertSizeRangeSD * getSD(datasetNumber))
			{
				// CP: the path length is within the insert size range
				listOfEdges.push_back(firstEdge);
				pathLengths.push_back(distanceOnLastEdge + pathLengths.at(level - 1));
				pathFound++;
				if(pathFound == 1)	// First path
				{

					for(UINT64 i = 0; i <listOfEdges.size() ; i++) // Copy the path.
						firstPath.push_back(listOfEdges.at(i));

					for(UINT64 i = 0; i <listOfEdges.size() - 1 ; i++) // All adjacent edges in the path are supported.
						flags.push_back(1);
				}
				else		// Not the first path.
				{
					UINT64 i , j;
					for( i = 0; i< firstPath.size() - 1; i++)		// Compare to the first path.
					{
						for( j = 0; j < listOfEdges.size() - 1; j++)
						{
							if(firstPath.at(i) == listOfEdges.at(j) &&  firstPath.at(i+1) == listOfEdges.at(j+1)) // A pair of edges adjacent in the first path is also adjacent in this path. We keep the support.
								break;
						}
						if(j == listOfEdges.size() - 1)		// A pair of adjacent edge in the first path is not adjacent in this path. So this pair is not supported anymore. So we clear the flag.
							flags.at(i) = 0;
					}
				}
				/*for(UINT64 i = 0; i <firstPath.size() ; i++) // Copy the path.
				{
					if(i <listOfEdges.size() - 1 )
						cout<< " (" << firstPath.at(i)->getSourceRead()->getReadNumber() <<", " <<firstPath.at(i)->getDestinationRead()->getReadNumber()<<") " <<flags.at(i);
					else
						cout<< " (" << firstPath.at(i)->getSourceRead()->getReadNumber() <<", " <<firstPath.at(i)->getDestinationRead()->getReadNumber()<<") ";
				}
				cout << endl;*/

				// CP: when found a valid path, return 1 path found and exit the recursive call
				return 1;
			}
			else		// add the new edge in the path
			{
				// CP: this length of this path is not within the valid range of the insert size and keep going deeper
				listOfEdges.push_back(firstEdge);
				pathLengths.push_back( distanceOnFirstEdge + pathLengths.at(level - 1) );
			}
		}
		else	// add the new edge in the path.
		{
			// CP: this path has not reached the destination edge yet and keep going deeper
			listOfEdges.push_back(firstEdge);
			pathLengths.push_back( distanceOnFirstEdge + pathLengths.at(level - 1) );
		}
	}

	for(UINT64 i = 0 ; i < graph->at(firstEdge->getDestinationRead()->getReadNumber())->size(); i++ )	// go deepter in the graph to reach destination edge
	{
		// CP: going to each edge that the last read of the firstEdge is connected to
		Edge * nextEdge =  graph->at(firstEdge->getDestinationRead()->getReadNumber())->at(i);
		// CP: if the two edges are compatible and the current path is not already too long, then call self again
		if(matchEdgeType(firstEdge, nextEdge) && pathLengths.at(level) < getMean(datasetNumber) + insertSizeRangeSD * getSD(datasetNumber))
			exploreGraph(nextEdge, lastEdge, nextEdge->getOverlapOffset(), distanceOnLastEdge, datasetNumber, level+1, firstPath, flags);
	}
	return pathFound;
}





/**********************************************************************************************************************
	Check for paths in the overlap graph between each matepair and calculate the support by using the paths.
**********************************************************************************************************************/
UINT64 OverlapGraph::findSupportByMatepairsAndMerge(void)
{
	CLOCKSTART;

	// if the file set is not mate-pair, then just skip
	if(meanOfInsertSizes.size() == 0) // no mate-pair
		return 0;
	// a path is defined by a vector of edges,
	vector <Edge *> copyOfPath;
	// This list of flags is used to mark the edges that are common in all paths found.
	vector <UINT64> copyOfFlags;
	UINT64 noPathsFound = 0, pathsFound = 0, mpOnSameEdge=0;

	vector <pairedEdges> listOfPairedEdges;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		Read * read1 = dataSet->getReadFromID(i);
		// initial mate-pair information is saved before building the graph using storeMatePairInfo()
		for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++)
		{
			Read * read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);
			// CP: read1 and read2 are paired-ends
			if(read1->getReadNumber() > read2->getReadNumber()) 		// To avoid finding path twice
				continue;
			// CP: ignore this pair if their dataset has an average insert size of 0
			if(meanOfInsertSizes.at(read1->getMatePairList()->at(j).datasetNumber) == 0)
				continue;

			if( findPathBetweenMatepairs(read1, read2, read1->getMatePairList()->at(j).matePairOrientation, read1->getMatePairList()->at(j).datasetNumber, copyOfPath, copyOfFlags) == true)
			{
				// Matepair on different edge
				if(copyOfPath.size() == 0)
					noPathsFound++;
				else
					pathsFound++;
			}
			else // Mate pair on the same edge
			{
				mpOnSameEdge++;
			}

			if(copyOfPath.size() > 1)	// Path found
			{
				for(UINT64 k = 0; k < copyOfFlags.size(); k++)
				{
					if(copyOfFlags.at(k) == 1)	// edge at k and k+1 is supported. We need to add it in the list if not present. If already the pair of edges present then increase the support
					{
						UINT64 l;
						for(l = 0; l < listOfPairedEdges.size(); l++)
						{
							if(listOfPairedEdges.at(l).edge1 == copyOfPath.at(k) && listOfPairedEdges.at(l).edge2 == copyOfPath.at(k+1)) // already present in the list
							{
								listOfPairedEdges.at(l).support = listOfPairedEdges.at(l).support + 1;	// only increase the support
								break;
							}
							else if(listOfPairedEdges.at(l).edge2->getReverseEdge() == copyOfPath.at(k) && listOfPairedEdges.at(l).edge1->getReverseEdge() == copyOfPath.at(k+1)) // already present in the list
							{
								listOfPairedEdges.at(l).support = listOfPairedEdges.at(l).support + 1;	// only increase the support
								break;
							}
						}
						if(l == listOfPairedEdges.size()) // not present in the list
						{
							if(copyOfPath.at(k)->getSourceRead()->getReadNumber() != copyOfPath.at(k)->getDestinationRead()->getReadNumber() || copyOfPath.at(k+1)->getSourceRead()->getReadNumber() != copyOfPath.at(k+1)->getDestinationRead()->getReadNumber()) // add in the list with support 1
							//if(copyOfPath.at(k)!=copyOfPath.at(k+1) && copyOfPath.at(k)!=copyOfPath.at(k+1)->getReverseEdge()) // do not want to add support between edge (a,a) and (a,a)
							{
								pairedEdges newPair;
								newPair.edge1 = copyOfPath.at(k);
								newPair.edge2 = copyOfPath.at(k+1);
								newPair.support = 1;
								newPair.isFreed = false;
								listOfPairedEdges.push_back(newPair);
							}
						}
					}
				}
			}
		}
	}

	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());

	UINT64 pairsOfEdgesMerged = 0;

	for(UINT64 i = 0; i<listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).support >= minimumSupport)// && listOfPairedEdges.at(i).edge1->flow > 0 && listOfPairedEdges.at(i).edge2->flow > 0)
		{
			pairsOfEdgesMerged++;
			cout << setw(4) << i + 1 << " Merging (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadNumber() << "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge1->flow << " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadNumber() << "," << setw(10) << listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge2->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->flow << " are supported " << setw(4) << listOfPairedEdges.at(i).support<<" times"<< endl;

			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();

			mergeEdges(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2);

			// BH: once we merge edge1 and edge2. We make sure that we do not try to merge these edges again. We mark all the pair of edes that contains edge1 and edge2 or their reverse edges.
			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r || listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r || listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}
	cout << pairsOfEdgesMerged <<" Pairs of Edges merged out of " << listOfPairedEdges.size() << " supported pairs of edges" <<endl;
	cout << "No paths found between " << noPathsFound << " matepairs that are on different edge." << endl;
	cout << "Paths found between " << pathsFound << " matepairs that are on different edge." << endl;
	cout << "Total matepairs on different edges " << pathsFound+ noPathsFound << endl;
	cout << "Total matepairs on same edge " << mpOnSameEdge << endl;
	cout << "Total matepairs " << pathsFound+noPathsFound+mpOnSameEdge << endl;
	CLOCKSTOP;
	return pairsOfEdgesMerged;

}


/**********************************************************************************************************************
	This function returns the string by overlapping the reads in an edge in the overlap graph
**********************************************************************************************************************/

string OverlapGraph::getStringInEdge(const Edge *edge)
{
	string read1, read2, readTemp, returnString;
	read1 = (edge->getOrientation() == 2 || edge->getOrientation() == 3) ?  edge->getSourceRead()->getStringForward() : edge->getSourceRead()->getStringReverse();
	read2 = (edge->getOrientation() == 1 || edge->getOrientation() == 3) ?  edge->getDestinationRead()->getStringForward() : edge->getDestinationRead()->getStringReverse();
	returnString = read1;
	UINT64 previousLength = read1.length(), substringLength;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++)
	{
		readTemp = (edge->getListOfOrientations()->at(i) == 1) ? dataSet->getReadFromID(edge->getListOfReads()->at(i))->getStringForward(): dataSet->getReadFromID(edge->getListOfReads()->at(i))->getStringReverse();

		substringLength =  readTemp.length() + edge->getListOfOverlapOffsets()->at(i) - previousLength;
		if( edge->getListOfOverlapOffsets()->at(i) ==  previousLength)
			returnString = returnString + 'N' + readTemp.substr(readTemp.length() - substringLength, substringLength);
		else
			returnString = returnString + readTemp.substr(readTemp.length() - substringLength, substringLength);
		previousLength = readTemp.length();


	}

	if(edge->getListOfReads()->empty()) // Simple edge
	{
		substringLength =  read2.length() + edge->getOverlapOffset() - read1.length();
		returnString = returnString + read2.substr(read2.length() - substringLength, substringLength);
	}
	else
	{
		substringLength = edge->getReverseEdge()->getListOfOverlapOffsets()->at(0);
		returnString = returnString + read2.substr(read2.length() - substringLength, substringLength);
	}
	return returnString;
}


/**********************************************************************************************************************
	This function removes in-trees and out-trees.
**********************************************************************************************************************/

UINT64 OverlapGraph::reduceTrees(void)
{
	CLOCKSTART;
	UINT64 NumOfInEdges, NumOfOutEdges, inFlow, outFlow, nodeMerged = 0;
	vector <Edge *> listOfInEdges, listOfOutEdges;
	for(UINT64 i = 0; i < graph->size(); i++)					// For each node in the graph
	{

		NumOfInEdges = 0; NumOfOutEdges = 0; inFlow = 0; outFlow = 0;
		listOfInEdges.clear(); listOfOutEdges.clear();
		for(UINT64 j = 0; j< graph->at(i)->size(); j++)
		{
			if(graph->at(i)->at(j)->flow == 0 || graph->at(i)->at(j)->flow != graph->at(i)->at(j)->getReverseEdge()->flow || graph->at(i)->at(j)->getSourceRead()->getReadNumber() == graph->at(i)->at(j)->getDestinationRead()->getReadNumber() )  // Some conditions for not being considered as a tree
					break;
			if(graph->at(i)->at(j)->getOrientation() == 0 || graph->at(i)->at(j)->getOrientation() == 1 )		// It is an in-edge
			{
				NumOfInEdges++;
				inFlow += graph->at(i)->at(j)->flow;															// Count the in flow
				listOfInEdges.push_back(graph->at(i)->at(j));
			}
			else																								// It is an out-edge
			{
				NumOfOutEdges++;
				outFlow += graph->at(i)->at(j)->flow;															// Count the out flow
				listOfOutEdges.push_back(graph->at(i)->at(j));
			}

			if(inFlow == outFlow && ( (NumOfInEdges == 1 && NumOfOutEdges > 1) || (NumOfInEdges > 1 && NumOfOutEdges == 1) ) )		// Either an in tree or an out tree
			{
				nodeMerged++;
				for(UINT64 k = 0; k < listOfInEdges.size(); k++)
				{
					for(UINT64 l = 0; l < listOfOutEdges.size(); l++)
					{
						mergeEdges(listOfInEdges.at(k)->getReverseEdge(), listOfOutEdges.at(l));
					}
				}
			}
		}
	}
	cout << setw(10) << nodeMerged << " trees removed." << endl;
	CLOCKSTOP;
	return nodeMerged;
}



/**********************************************************************************************************************
	removes composite path, dead-ends and trees
**********************************************************************************************************************/
bool OverlapGraph::simplifyGraph(void)
{
	UINT64 counter = 0;
	do
	{
		 counter  = removeDeadEndNodes();				// Remove dead-ends
		 counter += contractCompositePaths();		// Contract composite paths
		 counter += removeSimilarEdges();
		 counter += reduceTrees();					// Remove trees.
		 counter += reduceLoops();

	} while (counter > 0);
	return true;
}





/**********************************************************************************************************************
	Original Scaffolder function.
	*******************************************************************************************************************/
UINT64 OverlapGraph::scaffolder(void)
{
	CLOCKSTART;
	UINT64 pairsOfEdgesMerged = 0;
	vector <pairedEdges> listOfPairedEdges;

	// listOfCompositeEdges contains all the composite edges
	vector<Edge *> *listOfCompositeEdges = new vector<Edge *>;
	for(UINT64 i = 1; i < graph->size(); i++) // For each node
	{
		for(UINT64 j = 0 ; j < graph->at(i)->size(); j++) // for each edge
		{
			Edge *edge = graph->at(i)->at(j);
			if(!edge->getListOfReads()->empty()) // composite edge
			{
				listOfCompositeEdges->push_back(edge);	// List of all the composite edges.
			}
		}
	}

	for(UINT64 i = 0 ; i < listOfCompositeEdges->size(); i++) // For each composite edge in the graph
	{
		Edge *edge1 = listOfCompositeEdges->at(i);
		vector<Edge *> *listOfFeasibleEdges = getListOfFeasibleEdges(listOfCompositeEdges->at(i)); // Find the list of other edges that share unique matepairs with the current edge. Only check one of the endpoints of the current edge.
		for(UINT64 j = 0; j < listOfFeasibleEdges->size(); j++ ) // Check the current edge vs the list of edges for suppor for scaffolder
		{
			Edge *edge2 =listOfFeasibleEdges->at(j);
			UINT64 distance;
			UINT64 support = checkForScaffold(edge1,edge2,&distance); // check the support and distance
			if(support>0)	// If there are support then add the current pair in the list.
			{
				pairedEdges newPair;
				newPair.edge1 = edge1;
				newPair.edge2 = edge2;
				newPair.support = support;
				newPair.distance = distance;
				newPair.isFreed = false;
				listOfPairedEdges.push_back(newPair);

			}
		}
		delete listOfFeasibleEdges;	// Free the memory
	}
	delete listOfCompositeEdges;

	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());		// Sort the list according to support.

	for(UINT64 i = 0; i < listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).support >= minimumSupport)
		{
			pairsOfEdgesMerged++;
			listOfPairedEdges.at(i).distance = listOfPairedEdges.at(i).distance / listOfPairedEdges.at(i).support;
			cout << setw(4) << i + 1 << " (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadNumber() << "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge1->flow << " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadNumber() << "," << setw(10) << listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge2->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->flow << " are supported " << setw(4) << listOfPairedEdges.at(i).support << " times. Average distance: "<< setw(4) << listOfPairedEdges.at(i).distance << endl;
			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();
			mergeEdgesDisconnected(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2,listOfPairedEdges.at(i).distance);		// Merge the edges.
			// BH: if an edge is merged already, I make sure that I will not try to merge it again with other edges.
			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r || listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r || listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}

	CLOCKSTOP;
	return pairsOfEdgesMerged;
}

/**********************************************************************************************************************
 	 This functions returns a list of edges that might be joined to "edge"
 	 CP: For the input edge, find all the feasible edges that are linked with the input edge by a pair of edge-unique reads with appropriate distance
 	 CP: edge-unique reads are reads that are only present on one edge
***********************************************************************************************************************/
vector<Edge *> * OverlapGraph::getListOfFeasibleEdges(const Edge *edge)
{

	Edge * rEdge=edge->getReverseEdge(); // We want to find if there are other edges that share matepairs. current edge (u,v) we check the matepairs near the node v. That's why we took the reverse edge.
										// CP: why do you only check near v, not u?
										// BH: Here we are checkin if we can merge (u,v) followed by another edge. That why we only take the reads near v.
										// BH: At some point we will call this function with the reverse edge (v,u) as well. In that case we will look at the reads near vertex u.
	vector<Edge *> * feasibleListOfEdges = new vector<Edge *>;
	UINT64 dist = 0;
	for(UINT64 i = 0; i <rEdge->getListOfReads()->size(); i++) // for each read in the edge
	{
		// CP: dist is the distance from the beginning of rEdge to the beginning of this read.
		dist+=rEdge->getListOfOverlapOffsets()->at(i);	// offset. We do not have to go much deeper. we need to make sure that we atleast go upto the longest insert size.
		// CP: if this read is too far inside rEdge, then no need to go further inside. Stop and return the list
		if(dist > 2*longestMeanOfInsertSize)	// longest insert size mean
			 break;
		// CP: retrieve the current read: r1
		UINT64 mp1=rEdge->getListOfReads()->at(i); // mate pair 1
		Read *r1 = dataSet->getReadFromID(mp1); // read1

		if(r1->getListOfEdgesForward()->size() == 1) // only present in this current edge
		{
			// CP: r1 is a unique read that is present in this current edge only. Ignore non-unique reads
			// CP: this for loop considers r1 and forward edges of r2: r1->.......r2->
			for(UINT64 j = 0; j < r1->getMatePairList()->size(); j++) // for each matepair of current read1
			{
				// CP: r2 is the paired read of r1
				UINT64 mp2 = r1->getMatePairList()->at(j).matePairID; // matepair 2
				Read* r2 = dataSet->getReadFromID(mp2); // read2
				vector<Edge *> *list = r2->getListOfEdgesForward(); // location of read2
				// CP: use read2 if it's on one and only one edge and it's not on the input forward/reverse edge and its distance is adequate
				if(list->empty() || list->size() > 1 || list->at(0) == edge || list->at(0) == edge->getReverseEdge() || r2->getLocationOnEdgeForward()->at(0) > 2*longestMeanOfInsertSize) // Must be present uniquly on the edge and withing the distance of longest insert size.
					continue;
				UINT64 k;
				for(k = 0; k<feasibleListOfEdges->size();k++)		// add in the list of feasible edges. This list is expected to be small.
				{
					if(feasibleListOfEdges->at(k) == list->at(0))	// already in the list.
						break;
				}
				if(k == feasibleListOfEdges->size())	// Not presnet in the list.
				{
					feasibleListOfEdges->push_back(list->at(0));	// insert the edge in the list.
				}
			}
			// CP: this for loop considers r1 and reverse edges of r2: r1->.......r2<-
			// CP: what if the same read are on a forward edge or another reverse edge? Is this read still unique?
			for(UINT64 j = 0; j < r1->getMatePairList()->size(); j++)		// Same thing we do for the revese edges.
			{
				UINT64 mp2 = r1->getMatePairList()->at(j).matePairID;
				Read* r2 = dataSet->getReadFromID(mp2);
				vector<Edge *> *list = r2->getListOfEdgesReverse();
				if(list->empty() || list->size() > 1 || list->at(0) == edge || list->at(0) == edge->getReverseEdge() || r2->getLocationOnEdgeReverse()->at(0) > 2*longestMeanOfInsertSize)
					continue;
				UINT64 k;
				for(k = 0; k<feasibleListOfEdges->size();k++)
				{
					if(feasibleListOfEdges->at(k) == list->at(0))
						break;
				}
				if(k == feasibleListOfEdges->size())
				{
					feasibleListOfEdges->push_back(list->at(0));
				}
			}
		}
	}

	return feasibleListOfEdges;	// list of edges that might be joined with the current edge for scaffolding
}


/**********************************************************************************************************************
 	 Check how many unique matepairs support this pair of edges.
 	 CP: *edge1 and *edge2 are NOT modified here
 	 CP: return the number of matepairs supporting this pair or edges.
 	 CP: *distance (the distance between the ends of the two edges) is returned by reference
 	 CP: *distance is actually the sum of the distances measured by all matepairs. this would be easier to understand if it's the average distance
***********************************************************************************************************************/

UINT64 OverlapGraph::checkForScaffold(const Edge *edge1, const Edge *edge2,UINT64 *distance)
{
	UINT64 support = 0,dist = 0;
	*distance = 0;
	Edge *rEdge1 = edge1->getReverseEdge();
	vector<Edge *> *listRead1, *listRead2;		//  This is the lists of edges that contain read1 and read2
	vector<UINT64> *locationOnEdgeRead1, *locationOnEdgeRead2;
	// CP: listOfReads contains all the reads in the end section of edge1
	vector<UINT64> listOfReads;
	for(UINT64 i = 0; i <rEdge1->getListOfReads()->size(); i++)
	{
		dist+=rEdge1->getListOfOverlapOffsets()->at(i);
		if(dist>2*longestMeanOfInsertSize)
			 break;
		listOfReads.push_back(rEdge1->getListOfReads()->at(i));
	}
	for(UINT64 i = 0; i < listOfReads.size(); i++)
	{
		Read *read1 = dataSet->getReadFromID(listOfReads.at(i));
		for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++)// For each matepair
		{
			Read *read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);	// Get the read object of the matepair.
			//if(read1->getReadNumber() > read2->getReadNumber()) // To avoid duplicate computation
				//	continue;
			UINT64 orient = read1->getMatePairList()->at(j).matePairOrientation;		// Get the matepair orientation
			UINT64 datasetNumber = read1->getMatePairList()->at(j).datasetNumber;		// Get the dataset number

			// 0 = 00 means the reverse of r1 and the reverse of r2 are matepairs.
			// 1 = 01 means the reverse of r1 and the forward of r2 are matepairs.
			// 2 = 10 means the forward of r1 and the reverse of r2 are matepairs.
			// 3 = 11 means the forward of r1 and the forward of r2 are matepairs.
			// To calculate distance of forward read, flip the read and get the location of the offset.
			listRead1 = (orient == 0 || orient == 1) ? read1->getListOfEdgesForward() : read1->getListOfEdgesReverse();
			locationOnEdgeRead1 = (orient == 0 || orient == 1) ? read1->getLocationOnEdgeForward() : read1->getLocationOnEdgeReverse();
			// To calculate distance of reverse read, flip the read and get the location of the offset.
			listRead2 = (orient == 0 || orient == 2) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse();
			locationOnEdgeRead2 = (orient == 0 || orient == 2) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse();
			// Only consider uniquely mapped reads and the distance is less than mean+3*SD
			if( listRead1->size() == 1 && listRead2->size() == 1 && listRead1->at(0) == edge1->getReverseEdge() && listRead2->at(0) == edge2 && locationOnEdgeRead1->at(0) + locationOnEdgeRead2->at(0) < (getMean(datasetNumber) + insertSizeRangeSD * getSD(datasetNumber)) )  // Both the reads are present on only on edge and the distance is less that mean+3*sd
			{
				dist = locationOnEdgeRead1->at(0) + locationOnEdgeRead2->at(0);
				// if there are already in the same edge, don't do anything
				if(listRead1->at(0) == listRead2->at(0) ||  listRead1->at(0) == listRead2->at(0)->getReverseEdge()) // Not on the same edge
					continue;
				*distance +=dist;
				support++;
			}
		}
	}
	return support;
}





/**********************************************************************************************************************
	Original Scaffolder function.
	*******************************************************************************************************************/
UINT64 OverlapGraph::scaffolderTemp(void)
{
	CLOCKSTART;
	UINT64 pairsOfEdgesMerged = 0, dist;
	UINT8 orient, datasetNumber;
	Read *read1, *read2;
	vector <pairedEdges> listOfPairedEdges;
	vector<Edge *> *listRead1, *listRead2;
	vector<UINT64> *locationOnEdgeRead1, *locationOnEdgeRead2;

	for(UINT64 i = 1; i < graph->size(); i++)						// For each read
	{
		read1 = dataSet->getReadFromID(i);							// Get the read object
		for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++)// For each matepair
		{
			read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);	// Get the read object of the matepair.

			if(read1->getReadNumber() > read2->getReadNumber()) // To avoid duplicate computation
				continue;

			orient = read1->getMatePairList()->at(j).matePairOrientation;		// Get the matepair orientation
			datasetNumber = read1->getMatePairList()->at(j).datasetNumber;		// Get the dataset number

			// 0 = 00 means the reverse of r1 and the reverse of r2 are matepairs.
			// 1 = 01 means the reverse of r1 and the forward of r2 are matepairs.
			// 2 = 10 means the forward of r1 and the reverse of r2 are matepairs.
			// 3 = 11 means the forward of r1 and the forward of r2 are matepairs.
			// To calculate distance of forward read, flip the read and get the location of the offset.

			listRead1 = (orient == 0 || orient == 1) ? read1->getListOfEdgesForward() : read1->getListOfEdgesReverse();
			locationOnEdgeRead1 = (orient == 0 || orient == 1) ? read1->getLocationOnEdgeForward() : read1->getLocationOnEdgeReverse();

			// To calculate distance of reverse read, flip the read and get the location of the offset.
			listRead2 = (orient == 0 || orient == 2) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse();
			locationOnEdgeRead2 = (orient == 0 || orient == 2) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse();

			// Only consider uniquely mapped reads and the distance is less than mean+3*SD
			if( listRead1->size() == 1 && listRead2->size() == 1 && locationOnEdgeRead1->at(0) + locationOnEdgeRead2->at(0) < (getMean(datasetNumber) + insertSizeRangeSD * getSD(datasetNumber)) )  // Both the reads are present on only on edge and the distance is less that mean+3*sd
			{
				dist = locationOnEdgeRead1->at(0) + locationOnEdgeRead2->at(0);
				// if there are already in the same edge, don't do anything
				if(listRead1->at(0) == listRead2->at(0) ||  listRead1->at(0) == listRead2->at(0)->getReverseEdge()) // Not on the same edge
					continue;

				UINT64 k;

				for(k = 0; k < listOfPairedEdges.size(); k++)
				{
					if(listOfPairedEdges.at(k).edge1 == listRead1->at(0)->getReverseEdge() && listOfPairedEdges.at(k).edge2 == listRead2->at(0))		// This pair of edge was supported previously by another mate-pair. Only increase the support.
					{
						listOfPairedEdges.at(k).support += 1;
						listOfPairedEdges.at(k).distance += dist;
						break;
					}
					if(listOfPairedEdges.at(k).edge1 == listRead2->at(0)->getReverseEdge() && listOfPairedEdges.at(k).edge2 == listRead1->at(0))			// This pair of edge was supported previously by another mate-pair. Only increase the support.
					{
						listOfPairedEdges.at(k).support += 1;
						listOfPairedEdges.at(k).distance += dist;
						break;
					}
				}

				if(k == listOfPairedEdges.size())				// This mate pair is the first reads to support the pair of edges. Add the new pair of edges in the list with support 1.
				{
					pairedEdges newPair;
					newPair.edge1 = listRead1->at(0)->getReverseEdge();
					newPair.edge2 = listRead2->at(0);
					newPair.support = 1;
					newPair.distance = dist;
					newPair.isFreed = false;
					listOfPairedEdges.push_back(newPair);
				}
			}

		}
	}

	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());		// Sort the list according to support.

	for(UINT64 i = 0; i < listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).support >= minimumSupport)
		{
			pairsOfEdgesMerged++;
			listOfPairedEdges.at(i).distance = listOfPairedEdges.at(i).distance / listOfPairedEdges.at(i).support;
			cout << setw(4) << i + 1 << " (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadNumber() << "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge1->flow << " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadNumber() << "," << setw(10) << listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge2->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->flow << " are supported " << setw(4) << listOfPairedEdges.at(i).support << " times. Average distance: "<< setw(4) << listOfPairedEdges.at(i).distance << endl;
			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();

			mergeEdgesDisconnected(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2,listOfPairedEdges.at(i).distance);		// Merge the edges.

			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r || listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r || listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}

	CLOCKSTOP;
	return pairsOfEdgesMerged;
}

/**********************************************************************************************************************
	New Scaffolder function.
**********************************************************************************************************************/
/*
struct keyInfo
{
    Edge *edge1;   // key1 -> *edge1
    Edge *edge2;   // key2 -> *edge2

    keyInfo(Edge A, Edge B) :
        edge1(&A), edge2(&B) {}

    bool operator<(const keyInfo& A) const
    {
        return ((edge1<A.edge1) || (edge2<A.edge2)); 
    }

};

struct valueInfo
{
    UINT64 support;
    UINT64 distance;
    bool isFreed;

    valueInfo(UINT64 A,UINT64 B,bool C) : 
    support(A),distance(B),isFreed(C) {}
};

typedef map<keyInfo, valueInfo> MapType;


UINT64 OverlapGraph::scaffolder(void)
{
	CLOCKSTART;
	UINT64 pairsOfEdgesMerged = 0, dist;
	UINT8 orient, datasetNumber;
	Read *read1, *read2;
	vector <pairedEdges> listOfPairedEdges;
    // new map for struct key and val instead of vector
	MapType mapOfPairedEdges;
	vector<Edge *> *listRead1, *listRead2;
	vector<UINT64> *locationOnEdgeRead1, *locationOnEdgeRead2;

	for(UINT64 i = 1; i < graph->size(); i++)						// For each read
	{
		read1 = dataSet->getReadFromID(i);							// Get the read object
		for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++)// For each matepair
		{
			read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);	// Get the read object of the matepair.

			if(read1->getReadNumber() > read2->getReadNumber()) // To avoid duplicate computation
				continue;

			orient = read1->getMatePairList()->at(j).matePairOrientation;		// Get the matepair orientation
			datasetNumber = read1->getMatePairList()->at(j).datasetNumber;		// Get the dataset number

			// 0 = 00 means the reverse of r1 and the reverse of r2 are matepairs.
			// 1 = 01 means the reverse of r1 and the forward of r2 are matepairs.
			// 2 = 10 means the forward of r1 and the reverse of r2 are matepairs.
			// 3 = 11 means the forward of r1 and the forward of r2 are matepairs.
			// To calculate distance of forward read, flip the read and get the location of the offset.
			listRead1 = (orient == 0 || orient == 1) ? read1->getListOfEdgesForward() : read1->getListOfEdgesReverse();
			locationOnEdgeRead1 = (orient == 0 || orient == 1) ? read1->getLocationOnEdgeForward() : read1->getLocationOnEdgeReverse();
			// To calculate distance of reverse read, flip the read and get the location of the offset.
			listRead2 = (orient == 0 || orient == 2) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse();
			locationOnEdgeRead2 = (orient == 0 || orient == 2) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse();

			// get distance
			dist = locationOnEdgeRead1->at(0) + locationOnEdgeRead2->at(0);

			// Only consider uniquely mapped reads and the distance is less than mean+3*SD
			if( listRead1->size() == 1 && listRead2->size() == 1 && dist < (getMean(datasetNumber) + insertSizeRangeSD * getSD(datasetNumber)) )  // Both the reads are present on only on edge and the distance is less that mean+3*sd
			{
				// if there are already in the same edge, don't do anything
				if(listRead1->at(0) == listRead2->at(0) ||  listRead1->at(0) == listRead2->at(0)->getReverseEdge()) // Not on the same edge
					continue;

                // get key set (*edge1, *edge2)
                keyInfo keySet1(*listRead1->at(0)->getReverseEdge(), *listRead2->at(0));
                keyInfo keySet2(*listRead2->at(0)->getReverseEdge(), *listRead1->at(0));

                // find the keySet1 and keySet2
                MapType::iterator iter1=mapOfPairedEdges.find(keySet1);
                MapType::iterator iter2=mapOfPairedEdges.find(keySet2);

                // check keySet1
                if (iter1 != mapOfPairedEdges.end())        // if the keySet1 exists
                {
                    iter1->second.support += 1;             // add +1 the support
                    iter1->second.distance += dist;         // add +dist the distance
                }
                else if (iter2 != mapOfPairedEdges.end())   // if the keySet2 exists
                {
                    iter2->second.support += 1;             // add +1 the support
                    iter2->second.distance += dist;         // add +dist the distance
                }
                else                                        // if the keySet1 and keySet2 do not exist
                {
                    valueInfo valueSet(1, dist, false);     // get valueSet (support, distance, isFreed)
                    mapOfPairedEdges.insert(MapType::value_type(keySet1,valueSet)); // insert valueSet to the map
                }
			}
		}
	}

    


	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());		// Sort the list according to support.

	for(UINT64 i = 0; i < listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).support >= minimumSupport)
		{
			pairsOfEdgesMerged++;
			listOfPairedEdges.at(i).distance = listOfPairedEdges.at(i).distance / listOfPairedEdges.at(i).support;
			cout << setw(4) << i + 1 << " (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadNumber() << "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge1->flow << " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadNumber() << "," << setw(10) << listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge2->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->flow << " are supported " << setw(4) << listOfPairedEdges.at(i).support << " times. Average distance: "<< setw(4) << listOfPairedEdges.at(i).distance << endl;
			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();

			mergeEdgesDisconnected(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2,listOfPairedEdges.at(i).distance);		// Merge the edges.

			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r || listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r || listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}


	CLOCKSTOP;
	return pairsOfEdgesMerged;
}
*/

/**********************************************************************************************************************
	Check if two strings overlap. At least 10 bp must overlap.
	CP: return the length of the overlap. 0 if no overlap
**********************************************************************************************************************/

UINT64 OverlapGraph::findOverlap(const string & string1, const string & string2)
{
	UINT64 minimum = min(string1.length(), string2.length());
	for(UINT64 i = minimum - 1; i >= 10; i--)
	{
		if(string1.substr(string1.length() - i, i) == string2.substr(0,i))
		{
			return i;
		}
	}
	return 0;
}



/**********************************************************************************************************************
	Merge two edges that do not share any node.
**********************************************************************************************************************/
bool OverlapGraph::mergeEdgesDisconnected(Edge *edge1, Edge *edge2, UINT64 gapLength)
{
	if(edge1->getDestinationRead()->getReadNumber() == edge2->getSourceRead()->getReadNumber() && matchEdgeType (edge1, edge2)) // If the two edges are already connected. A --->B and B---->C. They share a common node
	{
		mergeEdges(edge1,edge2); // Merge the edges.
		return true;
	}

	// A------>B and C------>D. We need to check if the nodes B and C overlaps or not
	string string1, string2;
	string1 = ( edge1->getOrientation() == 1 || edge1->getOrientation() == 3 ) ? edge1->getDestinationRead()->getStringForward() : edge1->getDestinationRead()->getStringReverse(); // We will check if the two nodes overlap or not
	string2 = ( edge2->getOrientation() == 2 || edge2->getOrientation() == 3 ) ? edge2->getSourceRead()->getStringForward() : edge2->getSourceRead()->getStringReverse();

	UINT64 overlapLength = findOverlap(string1,string2); // Find the overlap between B and C. If they do not overlap then the return will be zero. We check for at least 10 bp overlap

	UINT64 overlapOffset1, overlapOffset2;
	Edge *edgeForward = new Edge();		// Forward edge
	Read *read1 = edge1->getSourceRead(), *read2 = edge2->getDestinationRead(); // Get the reads.
	Edge *edgeReverse = new Edge();		// Reverse edge.
	UINT8 orientationForward = mergedEdgeOrientationDisconnected(edge1,edge2); // Orientation of the forward edge based on the orientations of edge1 and edge2
	UINT8 orientationReverse = twinEdgeOrientation(orientationForward);			// Orientation of the reverse edge.

	if(overlapLength == 0) // Strings in the read B and C do not overlap
	{
		// CP: do you insert Ns if they don't overlap? It's important to insert Ns
		// BH: We do not add N here. But when we generate the string from the edges, we insert N's there.
		overlapOffset1 = edge1->getDestinationRead()->getReadLength();	// In this case we concatenate the strings in the edges. So the offset is the length of the read B
		overlapOffset2 = edge2->getSourceRead()->getReadLength();		// This is the overlap offset of the reverse edge.
	}
	else	// Strings in the read B and C do overlap
	{
		overlapOffset1 = edge1->getDestinationRead()->getReadLength() - overlapLength; // Overlap offset of the forward edge is taken according to the length of the string in B
		overlapOffset2 = edge2->getSourceRead()->getReadLength() - overlapLength; // overlap offset of the reverse edge is taken according to the length of the string in C
	}


	// CP: merge the forward edge
	vector<UINT64> * listReadsForward = new vector<UINT64>;		// List of reads in the forward edge.
	vector<UINT16> * listOverlapOffsetsForward= new vector<UINT16>;	// List of overlap offsets in the reads of the forward edge.
	vector<UINT8> * listOrientationsForward = new vector<UINT8>;	// List of orientations of the reads of the forward edge. 1 means forward string of the reads, 0 means reverse string of the read
	mergeListDisconnected(edge1, edge2, overlapOffset1, gapLength, listReadsForward, listOverlapOffsetsForward, listOrientationsForward);	// Merge the list of reads, overlaps etc for the forward edge
	edgeForward->makeEdge(read1,read2,orientationForward, edge1->getOverlapOffset() + edge2->getOverlapOffset() + overlapOffset1, listReadsForward, listOverlapOffsetsForward, listOrientationsForward);

	// CP: merge the reverse edge
	vector<UINT64> * listReadsReverse = new vector<UINT64>;			// List of reads in the reverse edge.
	vector<UINT16> * listOverlapOffsetsReverse= new vector<UINT16>;	// List of overlap offsets in the reads of the reverse edge.
	vector<UINT8> * listOrientationsReverse = new vector<UINT8>;	// List of orientations of the reads of the reverse edge. 1 means forward string of the reads, 0 means reverse string of the read
	mergeListDisconnected(edge2->getReverseEdge(),edge1->getReverseEdge(), overlapOffset2, gapLength, listReadsReverse, listOverlapOffsetsReverse,listOrientationsReverse); // Merge the list of reads, overlaps etc for the reverse edge
	// BH: lengthReverseEdge is the overlap offset of the new edge.
	UINT64 lengthReverseEdge = edge1->getReverseEdge()->getOverlapOffset() + edge2->getReverseEdge()->getOverlapOffset() + overlapOffset2;
	edgeReverse->makeEdge(read2, read1, orientationReverse, lengthReverseEdge, listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);

//	edgeReverse->makeEdge(read2, read1, orientationReverse, edge1->getReverseEdge()->getOverlapOffset() + edge2->getReverseEdge()->getOverlapOffset() + overlapOffset2, listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);

	edgeForward->setReverseEdge(edgeReverse);	// set the pointer of reverse edge
	edgeReverse->setReverseEdge(edgeForward);	// set the pointer of reverse edge

	UINT16 flow = min(edge1->flow,edge2->flow);	// Take the minimum of the flow from the two original edges.
	UINT64 coverage = min(edge1->coverageDepth, edge2->coverageDepth);	// not used
	edgeForward->flow = flow;	// set the flow in the forward edge.
	edgeForward->coverageDepth = coverage;

	edgeReverse->flow = flow;	// set the flow in the reverse edge.
	edgeReverse->coverageDepth = coverage;

	//if(flowComputed && flow == 0 && edgeForward->getOverlapOffset() > 1000)
	//{
	//	cout << "Check for error inserting edge between " << edgeForward->getSourceRead()->getReadNumber() << " and " << edgeForward->getDestinationRead()->getReadNumber() << " Length: " << edgeForward->getOverlapOffset() << endl;
	//}
	insertEdge(edgeForward); // insert forward the edge in the graph.
	insertEdge(edgeReverse); // insert the reverse edge in the graph.

	edge1->flow = edge1->flow - flow;		// decrease the flow in the original edge.
	edge1->getReverseEdge()->flow = edge1->getReverseEdge()->flow - flow; // decrease the flow in the original edge.
	edge1->coverageDepth = edge1->coverageDepth - coverage;
	edge1->getReverseEdge()->coverageDepth = edge1->getReverseEdge()->coverageDepth - coverage;

	edge2->flow = edge2->flow - flow; // decrease the flow in the original edge.
	edge2->getReverseEdge()->flow = edge2->getReverseEdge()->flow - flow; // decrease the flow in the original edge.
	edge2->coverageDepth = edge2->coverageDepth - coverage;
	edge2->getReverseEdge()->coverageDepth = edge2->getReverseEdge()->coverageDepth - coverage;

	if(edge1->flow == 0 || flow == 0)	// Remove the edge1 if the flow is used.
		removeEdge(edge1);
	if(edge2->flow == 0 || flow == 0)	// Remove the edge2 if the flow is used.
		removeEdge(edge2);

	return true;
}



/**********************************************************************************************************************
	Merge the list of reads, list of overlap offsets and list of orientations of two edges.
	CP: input: edge1 and edge2, NOT modified. All these pointers need to be changed to const Edge *edge1
	CP: the overlapOffset is the overlapOffset of the new edge, which is from the beginning of source read of the old edge1 to the beginning of the destination read of the old edge2.
	gapLength is not used here, but it's pairedEdges's distance
	CP: return-by-pointer: *listReads,  *listOverlapOffsets, *listOrientations
**********************************************************************************************************************/

bool OverlapGraph::mergeListDisconnected(Edge *edge1, Edge *edge2, UINT64 overlapOffset, UINT64 gapLength, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations)
{
	UINT64 sum = 0;

	// CP: Add all the reads in the first edge, NOT including its source read
	for(UINT64 i = 0; i < edge1->getListOfOrientations()->size(); i++)	// Get the list from the first edge.
	{
		listReads->push_back(edge1->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge1->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge1->getListOfOrientations()->at(i));
		sum += edge1->getListOfOverlapOffsets()->at(i);
	}

	// CP: Add the destination read of the first edge and set its offset and orientation
	listReads->push_back(edge1->getDestinationRead()->getReadNumber());		// CP: Add the destination read of the first edge.
	listOverlapOffsets->push_back(edge1->getOverlapOffset() - sum);			// CP: calculate the overlapOffset of the destination read
	if(edge1->getOrientation() == 1 || edge1->getOrientation() == 3)		// CP: calculate the orientation of the destination read
		listOrientations->push_back(1);
	else
		listOrientations->push_back(0);

	// CP: Add the source read of the second edge.
	listReads->push_back(edge2->getSourceRead()->getReadNumber());		// Add the source read of the second edge.
	listOverlapOffsets->push_back(overlapOffset);						// the overlapOff inputted from the function parameter
	if(edge2->getOrientation() == 2 || edge2->getOrientation() == 3)
		listOrientations->push_back(1);
	else
		listOrientations->push_back(0);

	// CP: Add all the reads in the second edge, NOT including its destination read
	for(UINT64 i = 0; i < edge2->getListOfOrientations()->size(); i++)	// Get the list from the second edge.
	{
		listReads->push_back(edge2->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge2->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge2->getListOfOrientations()->at(i));
	}
	return true;
}




/**********************************************************************************************************************
	Orientation of the merged edge
**********************************************************************************************************************/
UINT8 OverlapGraph::mergedEdgeOrientationDisconnected(const Edge *edge1, const Edge *edge2)
{
	UINT8 or1 = edge1->getOrientation(), or2 = edge2->getOrientation(),returnValue;
	     if( (or1 == 1 || or1 == 0) && (or2 == 0 || or2 == 2) )		// <---------*  and *-----------<
		returnValue = 0;
	else if( (or1 == 1 || or1 == 0) && (or2 == 1 || or2 == 3) )		// <---------*  and *----------->
		returnValue = 1;
	else if( (or1 == 2 || or1 == 3) && (or2 == 0 || or2 == 2) )		// >---------*  and *-----------<
		returnValue = 2;
	else if( (or1 == 2 || or1 == 3) && (or2 == 1 || or2 == 3) )		// >---------*  and *----------->
		returnValue = 3;
	else
	{
		cout<<(int)or1<<" "<<(int)or2<<endl;
		MYEXIT("Unable to merge.")
	}
	return returnValue;
}



/**********************************************************************************************************************
	Remove edges with similar endpoint in the overlap graph
**********************************************************************************************************************/

UINT64 OverlapGraph::removeSimilarEdges(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	vector <Edge *> listOfEdges1, listOfEdges2;
	vector <UINT64> listOfEditDistance;
	for(UINT64 i = 1; i < graph->size(); i++)	// For each node.
	{
		if(!graph->at(i)->empty())		// The node has some edges in the graph.
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)	// For all the edges.
			{
				Edge * e1 = graph->at(i)->at(j);
				UINT64 source1 = e1->getSourceRead()->getReadNumber(), destination1 = e1->getDestinationRead()->getReadNumber();
				if( source1 < destination1)
				{
					for(UINT64 k = j + 1; k < graph->at(i)->size(); k++)
					{
						Edge * e2 = graph->at(i)->at(k);
						UINT64 source2 = e2->getSourceRead()->getReadNumber(), destination2 = e2->getDestinationRead()->getReadNumber();
						if( source1 == source2 && destination1 == destination2)		// Check both edges starts and ends at same end points
						{
							if( abs((int)(e1->getOverlapOffset() - e2->getOverlapOffset())) <  e2->getOverlapOffset()/20 ) // The lengths are more than 95% similar
							{
								string s1 = getStringInEdge(e1);		// Get the string on the first edge.
								string s2 = getStringInEdge(e2);		// Get the string on the second edge.
								UINT64 editDistance = calculateEditDistance(s1,s2);		// Calculate the edit distance between the strings.
								if( editDistance < min(e1->getOverlapOffset(), e2->getOverlapOffset())/20 )	// If the edit distance is less than 5% of the length of the shortest string.
								{
									UINT64 l;
									for(l= 0; l <  listOfEdges1.size(); l++)	// Already in the list.
									{
										if(listOfEdges2.at(l) ==  e2 || listOfEdges2.at(l) == e1) // check if the edges are already used.
											break;
									}
									if(l ==  listOfEdges1.size())				// Not in the list. Add it in the list.
									{
										listOfEdges1.push_back(e1);				// We will keep this edge.
										listOfEdges2.push_back(e2);				// This edge will be deleted and the flow will be moved to the first edge.
										listOfEditDistance.push_back(editDistance);	// Also store the edit distance.
									}
								}
							}
						}
					}
				}
			}
		}
	}
	cout << listOfEdges1.size()<< " edges to remove" << endl;
	for(UINT64 i = 0; i < listOfEdges1.size(); i++)
	{
		cout << setw(10) << ++ counter << " removing edge ("<< setw(10) << listOfEdges1.at(i)->getSourceRead()->getReadNumber()<<"," << setw(10) << listOfEdges1.at(i)->getDestinationRead()->getReadNumber()<<") Lengths : " << setw(10) << listOfEdges1.at(i)->getOverlapOffset() << " and " << setw(10) << listOfEdges2.at(i)->getOverlapOffset() << " Flows: " << setw(3) << listOfEdges1.at(i)->flow << " and " << setw(3) << listOfEdges2.at(i)->flow << " Edit Distance: " << setw(5) << listOfEditDistance.at(i) << " Reads: " << listOfEdges1.at(i)->getListOfReads()->size() << " and " << listOfEdges2.at(i)->getListOfReads()->size() << endl;
		listOfEdges1.at(i)->flow += listOfEdges2.at(i)->flow;			// Move the flow of the delete edge to this edge.
		listOfEdges1.at(i)->getReverseEdge()->flow += listOfEdges2.at(i)->getReverseEdge()->flow;	// Same for the reverse edge.
		removeEdge(listOfEdges2.at(i));	// Then remove the edge with similar string.
	}
	cout << counter << " edges removed." << endl;
	CLOCKSTOP;
	return counter;
}



/**********************************************************************************************************************
	Resolve ambiguous nodes according to coverage.
**********************************************************************************************************************/
UINT64 OverlapGraph::resolveNodes(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	vector <Edge *> listOfInEdges, listOfOutEdges;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		listOfInEdges.clear();
		listOfOutEdges.clear();
		if(graph->at(i)->size()==4)								// Nodes with 4 edges.
		{
			for(UINT64 j=0; j < graph->at(i)->size(); j++)
			{
				if(graph->at(i)->at(j)->getSourceRead() == graph->at(i)->at(j)->getDestinationRead())	// Will not consider the nodes that has loop
				{
					listOfInEdges.clear();
					listOfOutEdges.clear();
					break;
				}
				if(graph->at(i)->at(j)->getOrientation() == 0 || graph->at(i)->at(j)->getOrientation() == 1)	// In-edges.
				{
					listOfInEdges.push_back(graph->at(i)->at(j)->getReverseEdge());
				}
				else																							// Out-edges
				{
					listOfOutEdges.push_back(graph->at(i)->at(j));
				}
			}
			if(listOfInEdges.size() == 2 && listOfOutEdges.size() == 2)											// If two in-edges and two out-edges.
			{

				Edge *inEdge1, *inEdge2, *outEdge1, *outEdge2;
				// Calculate the mean and SD of coverage depth of the unique reads in these edges.
				// CP: why not call this function for every edge of the graph at the beginning?
				// BH: I do not need this function for all the edges. I only called this function when I need the covrage depth.
				getBaseByBaseCoverage(listOfInEdges.at(0));
				getBaseByBaseCoverage(listOfInEdges.at(1));
				getBaseByBaseCoverage(listOfOutEdges.at(0));
				getBaseByBaseCoverage(listOfOutEdges.at(1));
				if(listOfInEdges.at(0)->coverageDepth > listOfInEdges.at(1)->coverageDepth)	// Sort the in-Edges.
				{
					inEdge1 = listOfInEdges.at(0);
					inEdge2 = listOfInEdges.at(1);
				}
				else
				{
					inEdge1 = listOfInEdges.at(1);
					inEdge2 = listOfInEdges.at(0);
				}

				if(listOfOutEdges.at(0)->coverageDepth > listOfOutEdges.at(1)->coverageDepth)	// Sort the out-edges.
				{
					outEdge1 = listOfOutEdges.at(0);
					outEdge2 = listOfOutEdges.at(1);
				}
				else
				{
					outEdge1 = listOfOutEdges.at(1);
					outEdge2 = listOfOutEdges.at(0);
				}

				UINT64 flag1 = 0, flag2 = 0;
				if(isOverlappintInterval(inEdge1->coverageDepth, inEdge1->SD, outEdge1->coverageDepth, outEdge1->SD) && !isOverlappintInterval(inEdge1->coverageDepth, inEdge1->SD, outEdge2->coverageDepth, outEdge2->SD) && !isOverlappintInterval(inEdge2->coverageDepth, inEdge2->SD, outEdge1->coverageDepth, outEdge1->SD))
				{
					flag1 = 1;
				}
				if(isOverlappintInterval(inEdge2->coverageDepth, inEdge2->SD, outEdge2->coverageDepth, outEdge2->SD) && !isOverlappintInterval(inEdge2->coverageDepth, inEdge2->SD, outEdge1->coverageDepth, outEdge1->SD) && !isOverlappintInterval(inEdge1->coverageDepth, inEdge1->SD, outEdge2->coverageDepth, outEdge2->SD))
				{
					flag2 = 1;
				}
				if(flag1 == 1 )
				{
					counter++;
					cout << setw(10) << counter << " Merging edges (" << setw(10) << inEdge1->getSourceRead()->getReadNumber() << "," << setw(10) <<inEdge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(6) << inEdge1->getOverlapOffset() << " Flow: " << setw(3) << inEdge1->flow << " Coverage: " << setw(4) << inEdge1->coverageDepth << " SD: " << setw(3) << inEdge1->SD << " and (" << setw(10) << outEdge1->getSourceRead()->getReadNumber() << "," << setw(10) << outEdge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(6) << outEdge1->getOverlapOffset() << " Flow: " << setw(3) << outEdge1->flow << " Coverage: " << setw(4) << outEdge1->coverageDepth << " SD: " << setw(3) << outEdge1->SD << endl;
					mergeEdges(inEdge1,outEdge1);		// Merge the edges.
				}
				if(flag2 == 1)
				{
					counter++;
					cout << setw(10) << counter << " Merging edges (" << setw(10) << inEdge2->getSourceRead()->getReadNumber() << "," << setw(10) <<inEdge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(6) << inEdge2->getOverlapOffset() << " Flow: " << setw(3) << inEdge2->flow << " Coverage: " << setw(4) << inEdge2->coverageDepth << " SD: " << setw(3) << inEdge2->SD << " and (" << setw(10) << outEdge2->getSourceRead()->getReadNumber() << "," << setw(10) << outEdge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(6) << outEdge2->getOverlapOffset() << " Flow: " << setw(3) << outEdge2->flow << " Coverage: " << setw(4) << outEdge2->coverageDepth << " SD: " << setw(3) << outEdge2->SD << endl;
					mergeEdges(inEdge2,outEdge2);			// Merge the edges.
				}
			}
		}
	}
	cout<< counter << " edges merged." << endl;
	CLOCKSTOP;
	return counter;
}




struct stackElement
{
	Edge * edge;
	UINT64 distance;
};

struct overlappingReads
{
	Edge *edge;
	UINT64 readID;
	UINT64 distance;
};



/**********************************************************************************************************************
	Calculate the coverage depth of an edge for every basepair and then update the Mean and SD of coverage depth in
	the edge. Only consider reads that are unique to the edge.
	CP: Calculate the variable, coverageDepth and SD, the Edge class
**********************************************************************************************************************/
void OverlapGraph::getBaseByBaseCoverage(Edge *edge)
{
	vector<UINT64> * coverageBaseByBase = new vector<UINT64>;
	UINT64 length = edge->getOverlapOffset() + edge->getDestinationRead()->getReadLength();		// Array lenght same as the string length in the edge.
	for(UINT64 i = 0; i <=length; i++)
	{
		coverageBaseByBase->push_back(0);				// At first all the bases are covered 0 times.
	}
	UINT64 overlapOffset = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++)	// For each read in the edge.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));
		overlapOffset += edge->getListOfOverlapOffsets()->at(i);		// Where the current read starts in the string.
		UINT64 readLength = read->getReadLength();
		for(UINT64 j = overlapOffset; j < overlapOffset + readLength; j++)
		{
			coverageBaseByBase->at(j) += read->getFrequency();		// Increase the coverage of all bases by the frequency of the read.
		}
	}

	overlapOffset = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) // Scan the reads again.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));
		overlapOffset += edge->getListOfOverlapOffsets()->at(i);
		UINT64 readLength = read->getReadLength();
		if(read->getListOfEdgesForward()->size() > 1)		// Clear the bases that are covered by reads apperaing in multiple places in the graph.
		{
			for(UINT64 j = overlapOffset; j < overlapOffset + readLength; j++)
			{
				coverageBaseByBase->at(j) = 0;
			}
		}
	}
	for(UINT64 i = 0; i < edge->getSourceRead()->getReadLength(); i++)		// For the source read, clear the bases. Because this read is present in multiple places or not all reads are considered for these bases.
	{
		coverageBaseByBase->at(i) = 0;
	}

	for(UINT64 i = 0; i < edge->getDestinationRead()->getReadLength(); i++)	// Similarly clear the bases covered by the destination read.
	{
		coverageBaseByBase->at(coverageBaseByBase->size() - 1 - i) = 0;
	}

	UINT64 sum = 0, variance=0, count = 0, mean = 0, sd = 0;

	for(UINT64 i = 0; i < coverageBaseByBase->size(); i++)
	{
		if(coverageBaseByBase->at(i) != 0)		// Count only the non-zero values.
		{
			sum += coverageBaseByBase->at(i);
			count++;
		}
	}
	if ( count != 0 )
	{
		mean = sum/count;		// Calculate the mean.

		for(UINT64 i = 0; i < coverageBaseByBase->size(); i++)
		{
			if(coverageBaseByBase->at(i) != 0)	// Calculate the variance.
			{
				variance += (mean - coverageBaseByBase->at(i)) * (mean - coverageBaseByBase->at(i));
			}
		}
		sd = sqrt(variance/count);	// Calculate the standard deviation.
	}
	edge->coverageDepth = mean;		// Update the mean of the current edge.
	edge->SD = sd;					// Update the standard deviation of the current edge.
	delete coverageBaseByBase;
}



/**********************************************************************************************************************
	For each node in the graph, sort all its incident edges according to destination read ID.
**********************************************************************************************************************/
void OverlapGraph::sortEdges()
{
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty())
		{
			sort(graph->at(i)->begin(), graph->at(i)->end(), compareEdgeID);
		}
	}
}

/**********************************************************************************************************************
	This function remove loops that can be traversed in only one way.
	a>--->b>--->b>--->c

	this function doesn't resolve loops that can be traversed in two ways:
	example: a<---<b<--->b>--->c
**********************************************************************************************************************/
UINT64 OverlapGraph::reduceLoops(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	Edge *ab,*bb,*bc;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(graph->at(i)->size() == 4) // only four edges. The loop is counted twice.
		{
			UINT64 loopCount = 0, incomingEdgeCount = 0, outgoingEdgeCount = 0;
			for(UINT64 j = 0; j< graph->at(i)->size(); j++)
			{
				if(graph->at(i)->at(j)->getDestinationRead()->getReadNumber() == i) // This is a loop
				{
					loopCount++;
					bb = graph->at(i)->at(j);
				}
				else if(graph->at(i)->at(j)->getOrientation() == 0 || graph->at(i)->at(j)->getOrientation() == 1) // incoming edge
				{
					incomingEdgeCount++;
					ab = graph->at(i)->at(j)->getReverseEdge();
				}
				else if(graph->at(i)->at(j)->getOrientation() == 2 || graph->at(i)->at(j)->getOrientation() == 3) // outgoing edge
				{
					outgoingEdgeCount++;
					bc = graph->at(i)->at(j);
				}
			}
			if(loopCount==2 && incomingEdgeCount == 1 && outgoingEdgeCount == 1)  // two in the loop and one incoming and one outgoing
			{
				cout<<"Loop found at node: " << i <<  " loop edge length: " << bb->getOverlapOffset() << " flow: " << bb->flow << " Other edge lengths: " << ab->getOverlapOffset() << " and " << bc->getOverlapOffset() << endl;
				//cout << getStringInEdge(bb) << endl;
				if(bb->getOrientation() == 0)
				{
					counter++;
					mergeEdges(ab,bb->getReverseEdge());
				}
				else if(bb->getOrientation() == 3)
				{
					counter++;
					mergeEdges(ab,bb);
				}
				else			// Arrow in the edge is >---< or <---->. In this case it is not possible to decide which way to go.
				{
					cout << "Unable to reduce loop because of the edge type." << endl;
				}
			}
		}
	}
	cout <<" Loops removed: " << counter << endl; // Number of loop we were able to reduce
	CLOCKSTOP;
	return counter;
}
