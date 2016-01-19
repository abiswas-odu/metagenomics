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
	Another constructor
**********************************************************************************************************************/
Read::Read(const string & s)
{
	// Initialize the variables.
	readNumber = 0;
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
	setRead(s);
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Read::~Read(void)
{
	// delete all the pointers.
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
	setFrequency(1);												// Set the frequency to 1.
	read = s;														// Store the string.
	readReverse = reverseComplement(s);
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
string Read::reverseComplement(const string & read)
{
	UINT64 stringLength = read.length();
	string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)						// Then complement the string. Change A to T, C to G, G to C and T to A.
	{
		if( read[i] & 0X02 ) // C or G
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
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




