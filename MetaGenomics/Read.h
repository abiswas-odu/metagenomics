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
		UINT64 fileIndex; 					// The sequence number of the read in the files...
		UINT64 readHashOffset; 					// The offset number of the read in the hash data...
		dna_bitset *read; 						// Binary representation of the read.
		string reverseComplement() const;
	public:
		UINT64 superReadID;						// 0 = not a contained read
												// otherwise superReadID contains the ID of the uniqe super read.
		Read(void);								// Default constructor.
		Read(const string & s, UINT64 fIndx);					// Another constructor.
		~Read(void);							// Destructor.

		bool setRead(const string & s); 		// Set the read.
		bool setReadNumber(UINT64 id); 			// Set the read number.
		void setFileIndex(UINT64 id); 			// Set the file index number.
		void setReadHashOffset(UINT64 offset); 			// Set the hash data table offset.
		void freeBitSet();				 			// Set the dna_bitset read pointer and save memory

		//bool compareReadOverlap(UINT64 seq1Start, UINT64 seq1Len, Read * seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient);
		string getStringFwd(void) const {return read->toString();} 									// Get the forward string of the current read.
		string getStringRev(void) const {return reverseComplement();} 								// Get the reverse string of the current read.
		UINT16 getReadLen(void) const {return read->getLength();} 								// Get the length of the string in the current read.

		UINT64 getReadNumber(void) const {return readNumber;} 								// Get the read number of the current read.
		UINT64 getFileIndex(void) const {return fileIndex;} 								// Get the fileIndex number of the current read.
		UINT64 getReadHashOffset(void)const {return readHashOffset;}  			// Get the hash data table offset.
};

#endif /* READS_H_ */
