/*
 * QueryDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYDATASET_H_
#define QUERYDATASET_H_

#include "Config.h"
#include "QueryRead.h"

class QueryDataset {

	UINT64 numberOfReads;								// Number of total reads present in the dataset.
	UINT64 numberOfUniqueReads; 						// number of unique reads in the dataset.

	UINT64 shortestReadLength;
	UINT64 longestReadLength;


	bool duplicateFilter();
	void sortReads();
//	bool compareReads (QueryRead *read1, QueryRead *read2);//for sortReads
public:
	QueryDataset();
	virtual ~QueryDataset();
	vector<QueryRead *> queryReadList;

	static bool qualityFilter(string & sequence);

	bool buildDataset(const string & QueryFilename);
	UINT64 getNumberOfReads(); 						// Get the number of total reads in the database.
	UINT64 getNumberOfUniqueReads(); 				// Get the number of unique reads in the database.

	QueryRead * getReadFromString(const string & read); 		// Find a read in the database given the string. Uses binary search in the list of reads.
	QueryRead * getReadFromID(UINT64 ID); 					// Find a read in the database given the ID in constant time.

	//for debugging
	bool printDataset(int from, int to); 					// Print few the reads in the dataset. For debuggin only.
	void saveReads(string fileName); // Save all the sorted unique reads in a text file. Used for debugging.

};

#endif /* QUERYDATASET_H_ */