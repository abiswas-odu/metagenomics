/*
 * Dataset.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef DATASET_H_
#define DATASET_H_
#include "Common.h"
#include "Read.h"


/**********************************************************************************************************************
	Class to store the dataset
**********************************************************************************************************************/
class Dataset
{
	private:
		UINT64 numberOfReads;								// Number of total reads present in the dataset.
		UINT64 numberOfUniqueReads; 						// number of unique reads in the dataset.
		UINT64 minimumOverlapLength;						// Length of the shortest read in the dataset.
		vector<Read *> *reads; 								// List of reads in the dataset.
		string reverseComplement(const string & read); 		// Get the reverse complement of a string.
		bool testRead(const string & read); 				// Test if the read is good. Contains only {A,C,G,T} and does not contain more than 80% of same base.
															// The dataset contains only good quality reads.
		bool removeDupicateReads(void); 					// Remove duplicate reads from the dataset. Frequency is stored for each of the reads.
		bool storeMatePairInformation(string fileName, UINT64 minOverlap, UINT64 datasetNumber); 	// Read the input file again and store mate pair information.
		void sortReads(void);								// Ted: sort the reads lexicographically.

	public:
		vector<string> pairedEndDatasetFileNames;
		vector<string> singleEndDatasetFileNames;
		UINT64 shortestReadLength;
		UINT64 longestReadLength;

		Dataset(void); 										// Default constructor.
		Dataset(vector<string> pairedEndFileNames, vector<string> singleEndFileNames, UINT64 minOverlap);// anotherconstructor.
		~Dataset(void);										// Default destructor.
		bool readDataset(string fileName, UINT64 minOverlap, UINT64 datasetNumber); // Read the database from fasta/fastq file. Matepairs should be one after another in the file.
		UINT64 getNumberOfReads(void); 						// Get the number of total reads in the database.
		UINT64 getNumberOfUniqueReads(void); 				// Get the number of unique reads in the database.
		bool printDataset(void); 							// Print few the reads in the dataset. For debuggin only.
        // Ted: change it to constant pointer
		Read * getReadFromString(const string & read); 		// Find a read in the database given the string. Uses binary search in the list of reads.
		Read * getReadFromID(UINT64 ID); 					// Find a read in the database given the ID in constant time.
		void readMatePairsFromFile(void);					// Read the matePairs from file
		void saveReads(string fileName);					// Save all the sorted unique reads in a text file. Used for debugging.

};


#endif /* DATASET_H_ */
