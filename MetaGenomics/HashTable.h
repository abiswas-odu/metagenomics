/*
 * HashTable.h
 *
 *  Created on: Apr 22, 2013
 *      Author: b72, b8b
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_
#include "Common.h"
#include "Dataset.h"

#define BIG_CONSTANT(x) (x##LLU)

/* useful constants */
enum
{
	BASE_A = 0x0,	/* binary: 00 */
	BASE_C = 0x1,	/*'binary: 01 */
	BASE_G = 0x2,	/* binary: 10 */
	BASE_T = 0x3,	/* binary: 11 */
};
#define BASE_MASK 0x0000000000000003	/* binary: 11 */

/**********************************************************************************************************************
	Class to store hashtable.
**********************************************************************************************************************/
class HashTable{
	private:
		Dataset *dataSet;							// Pointer of the dataset. this is NOT modified here
		MPI_Win win;
		vector<UINT64> *memoryHashPartitions;
		vector<UINT64> *memoryDataPartitions;
		vector<UINT64> *memoryReadCount;
		UINT64 hashTableSize; 						// Ted: Size of the hash table. This is the prime number of mod operation.
		UINT64 hashDataTableSize; 					// Ted: Size of the hash data table. This is based on the number of reads.
		UINT64 *hashTable; 							// AB: List of hash offset of each hash key location
		UINT64 *hashData; 							// AB: List of hash keys of read number (ID) and orientation.
		UINT16 hashStringLength;					// Ted: Length of prefix and suffix of the reads to hash. This is equal to the minumum overlap length.
		mutable UINT64 numberOfHashCollision;		// Counter to count the number of hash collisions. For debugging only.
													// It's mutable such that it can be modified in the const member function, getListOfReads
		bool insertIntoTable(Read *read, string forwardRead, UINT64 *hashDataLengths, int myid);	// Insert a string in the hash table.
		bool hashReadLengths(string forwardRead, UINT64 *hashRecordCounts); 					// Ted: Hash prefix and suffix of the read and its reverse complement in the hash table. Turn over to the constant
		void setHashTableSize(UINT64 size); 		// Set the size of the hash table.
		void setHashTableDataSize(int myid);		// Set the size of the hash data table.
		UINT64 getHashIndex(const string & subString) const;	//Convert a set of data bytes from hash data table to DNA string
		string toString(UINT64 hashDataIndex,UINT64 stringLen) const;
		UINT64 getReadLength(UINT64 globalOffset, int myid) const; 		// Get the length of the string in the read at offset. NOT threadsafe
	public:
		HashTable(UINT64 parallelProcessPoolSize);							// Default constructor.
		~HashTable();								// Destructor.
		void createOffsetTable();
		void insertDataset(Dataset *d, UINT64 minOverlapLength,UINT64 parallelThreadPoolSize,int myid);	// Insert the dataset in the hash table.
		vector<UINT64*> * setLocalHitList(const string readString, int myid); 			// Get the list of reads that contain subString as prefix or suffix.
		void deleteLocalHitList(vector<UINT64*> *localReadHits);
		map<UINT64,string> getLocalHitList(vector<UINT64*> *localReadHits, string subString, UINT64 subStringIndx) const;
		UINT64 hashFunction(const string & subString) const; 						// Hash function.
		UINT64 getHashTableSize(void) const {return hashTableSize;}		// Get the size of the hash table.
		UINT64 getHashStringLength() const {return hashStringLength;}		// Get the hash string length.
		Dataset * getDataset(void) const {return dataSet;}					// Get the pointer to the dataset.

		string getStringForward(UINT64 globalOffset, int myid) const; 			// Get the forward string of the read at offset.
		string getStringReverse(UINT64 globalOffset, int myid) const;  			// Get the reverse string of the read at offset.

		void readReadLengthsFromFile(string fileName, UINT64 minOverlap,UINT64 *hashRecordCounts);
		void populateReadLengths(UINT64 *hashRecordCounts);												//Populate the read lengths in the hash table for future offset calculation
		void populateReadData(int myid);												//Populate the read sequence in the hash data
		void readReadSequenceFromFile(string fileName, UINT64 minOverlap, UINT64 *hashDataLengths, UINT64 &readID, int myid);

		/*MPI Related Routines*/
		UINT64 getLocalOffset(UINT64 globalOffset, int myid) const;
		bool isGlobalOffsetInRange(UINT64 globalOffset, int myid) const;
		int getOffsetRank(UINT64 globalOffset) const;
		string toStringMPI(UINT64  *hashDataBlock,UINT64 stringLen,  UINT64 startOffset) const;
		UINT64 getMemoryReadCount(int myid);
		UINT64 getMaxMemoryReadCount();
		UINT64 getMemoryMaxLocalOffset(int rank);
		void endEpoch();
		UINT64 getLocalReadID(UINT64 localOffset, int myid) const;
		UINT64 getLocalReadLength(UINT64 localOffset, int myid) const;
		string getLocalStringForward(UINT64 localOffset, int myid) const;
		UINT64 getLocalNextOffset(UINT64 localOffset, int myid) const;
		UINT64 getLocalReadOrient(UINT64 localOffset, int myid) const;
		vector<UINT64*> * getReads(UINT64 startID, UINT64 endID, int myid);

		bool needsProcessing(UINT64 read1ID, string readString, int myid); 		// determine is the read mist be processed by this process or not return true of false

		void setLockAll(){ MPI_Win_lock_all(0, win); }
		void unLockAll(){ MPI_Win_unlock_all(win); }

};


#endif /* HASHTABLE_H_ */
