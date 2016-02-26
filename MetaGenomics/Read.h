/*
 * Reads.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef READ_H_
#define READ_H_

#include "Common.h"
#include "dna.h"

class Edge;

/**********************************************************************************************************************
	Class to store a read.
**********************************************************************************************************************/

class Read
{
	private:
		UINT64 readNumber; 						// Unique Identification of the read.
		string readName;
		dna_bitset *read; 							// String representation of the read.
		string reverseComplement() const;
	public:
		UINT64 superReadID;						// 0 = not a contained read
												// otherwise superReadID contains the ID of the uniqe super read.
		Read(void);								// Default constructor.
		Read(const string & s);					// Another constructor.
		~Read(void);							// Destructor.

		bool setRead(const string & s); 		// Set the read.
		bool setReadNumber(UINT64 id); 			// Set the read number.
		bool setReadName(string name);

		bool compareReadOverlap(UINT64 seq1Start, UINT64 seq1Len, Read * seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient);
		string getStringForward(void) const {return read->toString();} 									// Get the forward string of the current read.
		string getStringReverse(void) const {return reverseComplement();} 								// Get the reverse string of the current read.
		UINT16 getReadLength(void) const {return read->getLength();} 								// Get the length of the string in the current read.
		UINT64 getReadNumber(void) const {return readNumber;} 								// Get the read number of the current read.
		string getReadName(void) const {return readName;} 								// Get the read name of the current read.
		bool addMatePair(Read *r, UINT8 orientation, UINT64 datasetNumber);				// Add a matepair in the list.
};

#endif /* READS_H_ */
