/*
 * SingleKeyHashTable.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#include "SingleKeyHashTable.h"

SingleKeyHashTable::SingleKeyHashTable() {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength = Config::hashKeyLength;
	this->maxMismatch = Config::maxMismatch;
	this->numberOfMode = 4;
	this->dataSet = NULL;
	this->hashTable = NULL;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
}

SingleKeyHashTable::SingleKeyHashTable(QueryDataset * qDataset) {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength = Config::hashKeyLength;
	this->maxMismatch = Config::maxMismatch;
	this->numberOfMode = 4;
	this->dataSet = qDataset;
	this->hashTable = NULL;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
}

SingleKeyHashTable::~SingleKeyHashTable() {
	// TODO Auto-generated destructor stub
	if(this->hashTable!=NULL)
	delete hashTable;
	this->dataSet = NULL;//we don't delete dataSet here
}

bool SingleKeyHashTable::createHashTables()
{
	if(this->dataSet==NULL)
	{
		cout<<"no data set"<<endl;
		return false;
	}
	else if(this->hashTable!=NULL)
	{
		cout<<"Hash Table already exists."<<endl;
		return false;
	}
	else
	{
		this->hashTable = new HashTable(this->hashKeyLength,  this->dataSet, this->numberOfMode);

		return true;
	}
}

// 00 = 0 means prefix of the forward string.
// 01 = 1 means suffix of the forward string.
// 10 = 2 means prefix of the reverse string.
// 11 = 3 means suffix of the reverse string.
string SingleKeyHashTable::getReadSubstring(UINT64 readID, int mode)
{
	QueryRead * read = this->dataSet->getReadFromID(readID);
	string str = (mode == 0 || mode == 1) ? read->getSequence() : read->reverseComplement();
	string subStr = (mode == 0 || mode == 2) ? str.substr(0,this->hashKeyLength) : str.substr(str.length() - this->hashKeyLength, this->hashKeyLength);

	return subStr;
}


bool SingleKeyHashTable::insertQueryDataset(QueryDataset* d)
{
	if(this->hashTable==NULL)
	{
		cout<<"Hash Table hasn't been created yet."<<endl;
		return false;
	}
	else
	{
		UINT64 datasetsize = this->dataSet->getNumberOfUniqueReads();
//		omp_set_dynamic(0);
//		omp_set_num_threads(Config::numberOfThreads);
//		#pragma omp parallel
//			{
//		#pragma omp for schedule(dynamic)
		UINT64 currentID = 1;
		while(currentID<=datasetsize)
		{
			if(currentID%1000000 == 0)
				cout << setw(10) << currentID << " reads inserted in the hash table. " << endl;
			QueryRead * read = this->dataSet->getReadFromID(currentID);
			string forwardRead = read->getSequence();
			string reverseRead = read->reverseComplement();
			string prefixForward = forwardRead.substr(0,this->hashKeyLength); 											// Prefix of the forward string.
			string suffixForward = forwardRead.substr(forwardRead.length() - this->hashKeyLength,this->hashKeyLength);	// Suffix of the forward string.
			string prefixReverse = reverseRead.substr(0,this->hashKeyLength);											// Prefix of the reverse string.
			string suffixReverse = reverseRead.substr(reverseRead.length() - this->hashKeyLength,this->hashKeyLength);	// Suffix of the reverse string.
			insertQueryRead(read, prefixForward, 0);
			insertQueryRead(read, suffixForward, 1);
			insertQueryRead(read, prefixReverse, 2);
			insertQueryRead(read, suffixReverse, 3);
			currentID++;
		}
		cout<<"Hash Table "<<" maximum collision number is: "<< this->numberOfHashCollision<<endl;
		cout<<"Hash Table "<<" maximum single read collision number is: "<< this->maxSingleHashCollision<<endl;


//			}//end of parallel

		return true;

	}
}
bool SingleKeyHashTable::insertQueryRead(QueryRead *read, string subString, int mode)
{
	UINT64 currentCollision =0;

	UINT64 index = this->hashTable->hashFunction(subString);
	while(!hashTable->isEmptyAt(index))
	{
		map<int,vector<UINT64>*>::iterator p=hashTable->getDataVectorsAt(index)->begin();
		int keymode = p->first;
		UINT64 keyreadID = p->second->at(0);

		string keyStr = this->getReadSubstring(keyreadID,keymode);


		if(keyStr == subString)
				break;
		numberOfHashCollision++;
		currentCollision++;
		index = (index == hashTable->getHashTableSize() - 1) ? 0: index + 1; 	// Increment the index
	}
	hashTable->insertValueAt(index,mode,this->numberOfMode,read->getIdentifier());							// Add the string in the list.

	if(currentCollision> this->maxSingleHashCollision)
		this->maxSingleHashCollision = currentCollision;
	if(currentCollision > 1000)
	{
		cout << currentCollision << " collisions for read " << read->getIdentifier() << " " << subString << " " << mode << endl;
	}
	return true;
}

/*



bool SingleKeyHashTable::insertQueryDataset(QueryDataset* querydataset)
{

	queryDataSet = querydataset;
	UINT64 datasetsize = queryDataSet->getNumberOfUniqueReads();
	hashTableNameList.push_back("forwardprefix");
//	hashTableNameList.push_back("forwardsuffix");
//	hashTableNameList.push_back("reverseprefix");
//	hashTableNameList.push_back("reversesuffix");
	InitializeAllHashTables();

CLOCKSTART;
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(unsigned int i = 0; i< hashTableNameList.size(); i++)
		{
			string stringmode = hashTableNameList.at(i);
			UINT64 currentID = 1;
			while(currentID<=datasetsize)
			{
				if(currentID%1000000 == 0)
					cout << setw(10) << currentID << " reads inserted in the hash table. " << endl;
				QueryRead * queryread = queryDataSet->getReadFromID(currentID);
				insertQueryRead(queryread, stringmode);
				currentID++;
			}
			HashTable * currentHashTable = hashTableMap.at(stringmode);
			cout<<"Hash Table "<<stringmode<<" maximum collision number is: "<< currentHashTable->numberOfHashCollision<<endl;
			cout<<"Hash Table "<<stringmode<<" maximum single read collision number is: "<< currentHashTable->maxSingleHashCollision<<endl;
//--------
// clean up hash table
// ---------
			vector <DataVector *> *hashTableData = currentHashTable->hashTable;
			for(UINT64 k=0;k<hashTableData->size();k++)
			{
				DataVector *datavector = hashTableData->at(k);
				if(datavector->keystring=="")
				{
					delete datavector;
					hashTableData->at(k) = NULL;
				}
				else datavector->keystring = "";
			}
// ---------
		}
	}

CLOCKSTOP;
	return true;
}

bool SingleKeyHashTable::insertQueryRead(QueryRead *read, string mode)
{
	UINT64 readID = read->getIdentifier();
	string keystring = getReadSubstring(mode,readID);
	HashTable * currentHashTable = hashTableMap.at(mode);
	return  currentHashTable->insertIntoHashTable(keystring,readID);


}
*/

/*
bool SingleKeyHashTable::doAlignment(Alignment* align, string mode, int subjectStart)
{

	if(mode=="forwardprefix")
	{
		align->queryOrientation = true;
		if(align->subjectReadSequence.length()- subjectStart - hashKeyLength >= align->queryRead->getSequence().length() - hashKeyLength) // The overlap must continue till the end.
			return checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(subjectStart + hashKeyLength, align->subjectReadSequence.length()-(subjectStart + hashKeyLength));
		string restQuery = align->queryRead->getSequence().substr(hashKeyLength,  restSubject.length());
		int currentMismatchCount=0;
		for(unsigned int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
			}
		}
		align->subjectStart = -subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "forwardsuffix")
	{
		align->queryOrientation = true;
		if(align->queryRead->getSequence().length()-hashKeyLength <= subjectStart)
			return checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(0, subjectStart);
		string restQuery = align->queryRead->getSequence().substr(align->queryRead->getReadLength()-hashKeyLength-subjectStart, restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(align->queryRead->getReadLength()-hashKeyLength-subjectStart+i, restSubject.at(i)));
			}
		}
		align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "reverseprefix")
	{
		align->queryOrientation = false;
		if(align->subjectReadSequence.length()- subjectStart - hashKeyLength >= align->queryRead->getSequence().length() - hashKeyLength) // The overlap must continue till the end.
			return checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(subjectStart + hashKeyLength, align->subjectReadSequence.length()-(subjectStart + hashKeyLength));
		string restQuery = align->queryRead->reverseComplement().substr(hashKeyLength,  restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
			}
		}
		align->subjectStart = -subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "reversesuffix")
	{
		align->queryOrientation = false;
		if(align->queryRead->getSequence().length()-hashKeyLength <= subjectStart)
			return checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(0, subjectStart);
		string restQuery = align->queryRead->reverseComplement().substr(align->queryRead->getReadLength()-hashKeyLength-subjectStart, restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(align->queryRead->getReadLength()-hashKeyLength-subjectStart+i, restSubject.at(i)));
			}
		}
		align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else return false;

	return true;
}

//the choice of the start and stop position should meet the minimum overlap requirement.
bool SingleKeyHashTable::subjectWindowRange(int& startpoint, int& stoppoint, string mode, string& subjectRead)
{
	if(mode=="forwardprefix" || mode == "reverseprefix")
	{
		startpoint = 0;
		stoppoint = subjectRead.length()-this->minimumOverlapLength;
	}
	else if(mode == "forwardsuffix" || mode == "reversesuffix")
	{
		startpoint = this->minimumOverlapLength - hashKeyLength;
		stoppoint = subjectRead.length()-hashKeyLength;

	}
	else return false;

	return true;
}

bool SingleKeyHashTable::checkForContainedAlignment(Alignment* align, string mode, int subjectStart)
{
	string subjectString=align->subjectReadSequence; // Get the forward of read1
	string queryString="";
//	string queryString = (mode=="forwardprefix" || mode=="forwardsuffix") ? align->queryRead->getSequence() : align->queryRead->reverseComplement(); // Get the string in read2 based on the orientation.
	if(mode=="forwardprefix" || mode=="forwardsuffix")
	{
		queryString = align->queryRead->getSequence();
		align->queryOrientation = true;
	}
	else if(mode=="reverseprefix" || mode=="reversesuffix")
	{
		queryString = align->queryRead->reverseComplement();
		align->queryOrientation = false;
	}
	else return false;

	if(mode=="forwardprefix" || mode=="reverseprefix")
									// mode = forwardprefix
									//   >--------MMMMMMMMMMMMMMM*******------> subject read1      M means match found by hash table
									//            MMMMMMMMMMMMMMM*******>       query read2      * means we need to check these characters for match
									//				OR
									// mode = reverseprefix
									//	 >---*****MMMMMMMMMMMMMMM*******------> subject read1
									//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of query read2
	{
		int restSubjectLength = subjectString.length() - subjectStart - this->hashKeyLength; 	// This is the remaining of read1
		int restQueryLength = queryString.length() - this->hashKeyLength; 	// This is the remaining of read2
		if(restSubjectLength >= restQueryLength)
		{
			string restSubject = subjectString.substr(subjectStart + hashKeyLength, restQueryLength);
			string restQuery = queryString.substr(hashKeyLength,  restQueryLength);
			int currentMismatchCount=0;
			for(int i=0; i<restQuery.length();i++)
			{
				if(restQuery.at(i)!=restSubject.at(i))
				{
					currentMismatchCount++;
					if(currentMismatchCount>maxMismatch)return false;
					align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
				}
			}
			align->subjectStart = -subjectStart;
			align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
			align->queryEnd = align->queryRead->getReadLength()-1;

		}
	}
	else if(mode=="forwardsuffix" || mode=="reversesuffix")
									// mode = forwardsuffix
									//   >---*****MMMMMMMMMMMMMMM-------------> subject read1      M means match found by hash table
									//      >*****MMMMMMMMMMMMMMM       		query read2      * means we need to check these characters for match
									//				OR
									// mode = reversesuffix
									//	 >---*****MMMMMMMMMMMMMMM-------------> subject read1
									//		<*****MMMMMMMMMMMMMMM				Reverse Complement of query Read2
	{
		int restSubjectLength = subjectStart;
		int restQueryLength = queryString.length() - this->hashKeyLength;
		if(restSubjectLength >= restQueryLength)
		{
			string restSubject = subjectString.substr(subjectStart-restQueryLength, restQueryLength);
			string restQuery = queryString.substr(0, restQueryLength);
			int currentMismatchCount=0;
			for(int i=0; i<restQuery.length();i++)
			{
				if(restQuery.at(i)!=restSubject.at(i))
				{
					currentMismatchCount++;
					if(currentMismatchCount>maxMismatch)return false;
					align->editInfor.insert(std::pair<int, char>(i, restSubject.at(i)));
				}
			}
			align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
			align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
			align->queryEnd = align->queryRead->getReadLength()-1;
		}
	}
	else return false;

	return true;

}

bool SingleKeyHashTable::singleKeySearch(edge & Edge)
{

	string subjectRead = Edge.subjectReadSequence;


	for(int i =0; i<this->hashTableNameList.size();i++)
	{


	string modestring = this->hashTableNameList.at(i);

	int startpoint,stoppoint;
	if(!subjectWindowRange(startpoint, stoppoint, modestring, subjectRead)) return false;//guarantee it meets the minimum overlaplength requirement for alignments.

	for(int j=startpoint;j<=stoppoint;j++)
	{
	string subString = subjectRead.substr(j, hashKeyLength);
	vector<UINT64> * currentIDList = hashTableMap.at(modestring)->getReadIDListOfReads(subString);
	for(UINT64 k=0;currentIDList!=NULL&&k<currentIDList->size();k++)
	{
		UINT64 currentID = currentIDList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName>=Edge.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		//however, this will cause missing of the contained read detection/alignment. Accept this reality.
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = Edge.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(doAlignment(align, modestring, j)==false) delete align;
			else Edge.alignmentList.push_back(align);

		}
	}
	}

	}

}

bool SingleKeyHashTable::singleKeySearch(SubjectAlignment & subjectAlign)
{

	string subjectRead = subjectAlign.subjectReadSequence;


	for(int i =0; i<this->hashTableNameList.size();i++)
	{


	string modestring = this->hashTableNameList.at(i);

	int startpoint,stoppoint;
	if(!subjectWindowRange(startpoint, stoppoint, modestring, subjectRead)) return false;//guarantee it meets the minimum overlaplength requirement for alignments.

	for(int j=startpoint;j<=stoppoint;j++)
	{
	string subString = subjectRead.substr(j, hashKeyLength);
	vector<UINT64>* currentIDList = hashTableMap.at(modestring)->getReadIDListOfReads(subString);
	for(UINT64 k=0;currentIDList!=NULL&&k<currentIDList->size();k++)
	{
		UINT64 currentID = currentIDList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName==subjectAlign.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		//however, this will cause missing of the contained read detection/alignment. Accept this reality.
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = subjectAlign.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(doAlignment(align, modestring, j)==false) delete align;
			else subjectAlign.queryAlignmentList.push_back(align);

		}
	}
	}

	}

}
*/
