/*
 * Reads.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef READ_H_
#define READ_H_

#include "Common.h"

class Edge;

struct MPlist
{
	UINT64 matePairID;				// This list stores the list of ID's of the matepairs of the current read.
	UINT8 matePairOrientation; 		// List of mate pair orientation.
									// 0 = 00 means the reverse of this read and the reverse of the second read are matepairs.
									// 1 = 01 means the reverse of this read and the forward of the second read are matepairs.
									// 2 = 10 means the forward of this read and the reverse of the second read are matepairs.
									// 3 = 11 means the forward of this read and the forward of the second read are matepairs.
	UINT8 datasetNumber;			//Dataset number.
};

/**********************************************************************************************************************
	Class to store a read.
**********************************************************************************************************************/

class Read
{
	private:
		UINT64 readNumber; 						// Unique Identification of the read.
		string read; 							// String representation of the read.
		string readReverse;
		UINT32 frequency; 						// Frequency of the read. Number of times this read is present in the dataset. Used in some statistical analysis.
		vector<MPlist> *matePairList;
		vector<Edge *> *listOfEdgesForward;   	// List of edges that contain the forward string of this read.
		vector<UINT64> *locationOnEdgeForward;	// List of locations on the edges that contain the forward string of the current read.
		vector<Edge *> *listOfEdgesReverse;		// List of edges that contain the reverse string of this read.
		vector<UINT64> *locationOnEdgeReverse;	// List of locations on the edges that contain the reverse string of the current read.

		string reverseComplement(const string & read);

	public:
		UINT64 coverageDepth;					// Estimated depth of coverage.  --> Not use anymore
		UINT64 locationInDataset;				// Not for debugging only. Will be removed.  --> Not use anymore
		bool isContainedRead;
		UINT64 superReadID;						// 0 = not a contained read
												// otherwise superReadID contains the ID of the uniqe super read.
		Read(void);								// Default constructor.
		Read(const string & s);					// Another constructor.
		~Read(void);							// Destructor.

		bool setRead(const string & s); 		// Set the read.
		bool setReadNumber(UINT64 id); 			// Set the read number.
		bool setFrequency(UINT32 freq);			// Set the ferquency of the read.

		string getStringForward(void){return read;} 									// Get the forward string of the current read.
		string getStringReverse(void){return readReverse;} 								// Get the reverse string of the current read.
		UINT16 getReadLength(void){return read.length();} 								// Get the length of the string in the current read.
		UINT64 getReadNumber(void) {return readNumber;} 								// Get the read number of the current read.
		UINT32 getFrequency(void) {return frequency;}									// Get the frequency of the current read.
		vector<MPlist> * getMatePairList(void) {return matePairList;} 					// Get the list of matepairs.
		vector<Edge *> * getListOfEdgesForward(void){return listOfEdgesForward;}		// Get the list of edges that contain the forward string of the current read.
		vector<UINT64> * getLocationOnEdgeForward(void){return locationOnEdgeForward;}	// Get the list of locations on the edges that contain the forward string of the current read.
		vector<Edge *> * getListOfEdgesReverse(void){return listOfEdgesReverse;}		// Get the list of edges that contain the reverse string of the current read.
		vector<UINT64> * getLocationOnEdgeReverse(void){return locationOnEdgeReverse;}	// Get the list of locations on the edges that contain the reverse string of the current read.
		bool addMatePair(Read *r, UINT8 orientation, UINT64 datasetNumber);				// Add a matepair in the list.

};

#endif /* READS_H_ */
