/*
 * Read.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#include "Common.h"
#include "Read.h"



/**********************************************************************************************************************
	Default constructor
**********************************************************************************************************************/
Read::Read(void)
{
	// Initialize the variables.
	readNumber = 0;
	readName="";
	frequency = 0;
	isContainedRead = false;
	superReadID = 0;
	matePairList = new vector<MPlist>;
	matePairList->resize(matePairList->size());						// Resize to 0 to reduce space.
	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	listOfEdgesReverse = new vector<Edge *>;
	listOfEdgesReverse->resize(listOfEdgesReverse->size());			// Resize to 0 to reduce space.
	locationOnEdgeForward = new vector<UINT64>;
	locationOnEdgeForward->resize(locationOnEdgeForward->size());	// Resize to 0 to reduce space.
	locationOnEdgeReverse = new vector<UINT64>;
	locationOnEdgeReverse->resize(locationOnEdgeReverse->size());	// Resize to 0 to reduce space.
	read = new dna_bitset();
}

/**********************************************************************************************************************
	Another constructor
**********************************************************************************************************************/
Read::Read(const string & s)
{
	// Initialize the variables.
	read = new dna_bitset(s.c_str(), s.length());
	setFrequency(1);
	readNumber = 0;
	readName="";
	frequency = 0;
	isContainedRead = false;
	superReadID = 0;
	matePairList = new vector<MPlist>;
	matePairList->resize(matePairList->size());						// Resize to 0 to reduce space.
	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	listOfEdgesReverse = new vector<Edge *>;
	listOfEdgesReverse->resize(listOfEdgesReverse->size());			// Resize to 0 to reduce space.
	locationOnEdgeForward = new vector<UINT64>;
	locationOnEdgeForward->resize(locationOnEdgeForward->size());	// Resize to 0 to reduce space.
	locationOnEdgeReverse = new vector<UINT64>;
	locationOnEdgeReverse->resize(locationOnEdgeReverse->size());	// Resize to 0 to reduce space.
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Read::~Read(void)
{
	// delete all the pointers.
	delete read;
	delete matePairList;
	delete listOfEdgesForward;
	delete listOfEdgesReverse;
	delete locationOnEdgeForward;
	delete locationOnEdgeReverse;

}

/**********************************************************************************************************************
	Function to store the read.
**********************************************************************************************************************/
bool Read::setRead(const string & s)
{
	// Ted: if(s.length() < 10) MYEXIT("Length of string less than 10.");	// Reads must be at least 10 bases in length.  --> not need anymore
	setFrequency(1);									// Set the frequency to 1.
	dna_bitset *tmpRead = new dna_bitset(s.c_str(), s.length());
	delete read;
	read=NULL;
	read = tmpRead;
	return true;
}
/**********************************************************************************************************************
	This function assigns an ID to the read.
**********************************************************************************************************************/
bool Read::setReadNumber(UINT64 id)
{
	if(id <= 0) MYEXIT("ID less than 1.");
	readNumber = id;												// Set the read number.
	return true;
}

/**********************************************************************************************************************
	This function assigns an Name to the read.
**********************************************************************************************************************/
bool Read::setReadName(string name)
{
	readName = name;
	return true;
}


/**********************************************************************************************************************
	This function sets the frequency of the read.
**********************************************************************************************************************/
bool Read::setFrequency(UINT32 freq)
{
	if(freq < 1) MYEXIT("Frequency less than 1.");
	frequency = freq;												// Set the frequency of the read.
	return true;
}

/**********************************************************************************************************************
	Returns the reverse complement of a read.
**********************************************************************************************************************/
string Read::reverseComplement() const {
	return read->toRevComplement();
}
/**********************************************************************************************************************
	This function adds a matpair
**********************************************************************************************************************/
bool Read::addMatePair(Read *r, UINT8 orientation, UINT64 datasetNumber)
{
	UINT64 ID = r->getReadNumber(), i;
	if(matePairList->empty()) 							//matePairList is empty.
	{
		MPlist newMP;
		newMP.matePairID = ID;
		newMP.matePairOrientation = orientation;
		newMP.datasetNumber = datasetNumber;
		matePairList->push_back(newMP);
	}
	else	// matePairList is not empty.
	{
		for(i = 0; i < matePairList->size(); i++) // Check if the matepair already present in the list.
		{
			if(matePairList->at(i).matePairID == ID && matePairList->at(i).matePairOrientation == orientation && matePairList->at(i).datasetNumber == datasetNumber) // Already present in the list.
			{
				//cout<< this->getReadNumber()<< " " << r->getReadNumber() << " already present in the dataset" << endl;
				//cout<< this->getStringForward() << " " << r->getStringForward() << endl;
				break;
			}
		}
		if (i == matePairList->size()) // matePairList does not contain ID and orientation in the list
		{
			MPlist newMP;
			newMP.matePairID = ID;
			newMP.matePairOrientation = orientation;
			newMP.datasetNumber = datasetNumber;
			matePairList->push_back(newMP);
		}
	}
	matePairList->resize(matePairList->size());					// Resize the list to reduce space.

	return true;
}

bool Read::compareReadOverlap(UINT64 seq1Start, UINT64 seq1Len, Read * seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient)
{
	return read->compareSubString(seq1Start,seq1Len,seq2->read,seq2Start,seq2Len, orient);
}




