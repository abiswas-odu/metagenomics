/*
 * Dataset.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#include "Common.h"
#include "Dataset.h"



/**********************************************************************************************************************
	Function to compare two reads. Used for sorting.
**********************************************************************************************************************/
bool compareReads (Read *read1, Read *read2)
{
	return read1->getStringForward() < read2->getStringForward();
}

/**********************************************************************************************************************
	Default constructor
**********************************************************************************************************************/
Dataset::Dataset(void)
{
	// Initialize the variables.
	numberOfUniqueReads = 0;
	numberOfReads = 0;
	minimumOverlapLength = 0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
	reads = new vector<Read *>;

}


/**********************************************************************************************************************
	Another constructor
**********************************************************************************************************************/
Dataset::Dataset(vector<string> pairedEndFileNames, vector<string> singleEndFileNames, UINT64 minOverlap)
{
	// Initialize the variables.
	numberOfUniqueReads = 0;
	numberOfReads = 0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
	reads = new vector<Read *>;
	pairedEndDatasetFileNames = pairedEndFileNames;
	singleEndDatasetFileNames = singleEndFileNames;
	minimumOverlapLength = minOverlap;
	UINT64 counter = 0;

	for(UINT64 i = 0; i < pairedEndDatasetFileNames.size(); i++)						// Read the paired-end datasets.
	{
		readDataset(pairedEndDatasetFileNames.at(i), minimumOverlapLength, counter++);
	}

	for(UINT64 i = 0; i < singleEndDatasetFileNames.size(); i++)						// Read the single-end datasets.
	{
		readDataset(singleEndDatasetFileNames.at(i), minimumOverlapLength, counter++);
	}
	cout << endl << "Shortest read length in all datasets: " << setw(5) << shortestReadLength<<endl;
	cout << " Longest read length in all datasets: " << setw(5) << longestReadLength <<endl;

	for(UINT64 i = 0 ; i < reads->size(); i++) 		// Assing ID's to the reads.
		reads->at(i)->setReadNumber(i + 1);
	numberOfUniqueReads=reads->size();
	reads->resize(reads->size());

	//sortReads();
	//removeDupicateReads();								// Remove duplicated reads for the dataset.

}


/**********************************************************************************************************************
	This function saves the lexicographically sorted unique reads in a file. Used for debugging only.
**********************************************************************************************************************/
void Dataset::saveReads(string fileName)
{
	ofstream outputFile;
	outputFile.open(fileName.c_str());
	if(outputFile == NULL)
		MYEXIT("Unable to open file: "+fileName);
	for(UINT64 i = 1; i <= numberOfUniqueReads; i++)
	{
		Read * read1 = getReadFromID(i);
		if(read1->superReadID!=0)
		{
			outputFile << setw(10) << i << " Contained in " <<  setw(10) << read1->superReadID << " " << read1->getStringForward() << endl;
		}
		else
		{
			outputFile << setw(10) << i << " Noncontained " <<  setw(10) << read1->superReadID << " " << read1->getStringForward() << endl;
		}
	}
	outputFile.close();
}

/**********************************************************************************************************************
	This function reads the dataset from FASTA/FASTQ files
**********************************************************************************************************************/
bool Dataset::readDataset(string fileName, UINT64 minOverlap, UINT64 datasetNumber)
{
	CLOCKSTART;
	cout << "Reading dataset: " << datasetNumber << " from file: " << fileName << endl;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
		MYEXIT("Unable to open file: "+fileName)
	UINT64 goodReads = 0, badReads = 0;
	vector<string> line;
	string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
	enum FileType { FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	while(!myFile.eof())
	{
		if( (goodReads + badReads ) != 0 && (goodReads + badReads)%1000000 == 0)
			cout<< setw(10) << goodReads + badReads << " reads processed in dataset " << setw(2) << datasetNumber <<". " << setw(10) << goodReads << " good reads." << setw(10) << badReads << " bad reads." << endl;
		if(fileType == UNDEFINED)
		{
			getline (myFile,text);
			if(text[0] == '>')
				fileType = FASTA;
			else if(text[0] == '@')
				fileType = FASTQ;
			else
				MYEXIT("Unknown input file format.");
			myFile.seekg(0, ios::beg);
		}
		line.clear();
		if(fileType == FASTA) 			// Fasta file
		{
			getline (myFile,text);
			line.push_back(text);
			getline (myFile,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line.at(0).erase(std::remove(line.at(0).begin(),line.at(0).end(),'>'),line.at(0).end());			//Sequence name
			line0=line.at(0);
			line1 = line.at(1);								// The first string is in the 2nd line.

		}
		else if(fileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (myFile,text);
				line.push_back(text);
			}
			line.at(0).erase(std::remove(line.at(0).begin(),line.at(0).end(),'>'),line.at(0).end());			//Sequence name
			line0=line.at(0);
			line1 = line.at(1); 			// The first string is in the 2nd line.
		}
		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > minOverlap && testRead(line1) ) // Test the read is of good quality.
		{

			line1ReverseComplement = reverseComplement(line1);
			if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
				line1ReverseComplement= line1;
			Read *r1=new Read(line1ReverseComplement);
			UINT64 len = r1->getReadLength();
			r1->setReadName(line0);
			if(len > longestReadLength)
				longestReadLength = len;
			if(len < shortestReadLength)
				shortestReadLength = len;

			reads->push_back(r1);						// Store the first string in the dataset.
			numberOfReads++;							// Counter of the total number of reads.
			goodReads++;
		}
		else
			badReads++;
	}

	myFile.close();
    cout << endl << "Dataset: " << setw(2) << datasetNumber << endl;
    cout << "File name: " << fileName << endl;
	cout << setw(10) << goodReads << " good reads in current dataset."  << endl;
	cout << setw(10) << badReads << " bad reads in current dataset." << endl;
	cout << setw(10) << goodReads + badReads << " total reads in current dataset." << endl;
	cout << setw(10) << numberOfReads << " good reads in all datasets." << endl << endl;;
	CLOCKSTOP;
	return true;
}


// Ted: multi-thread this sort later
void Dataset::sortReads(void)
{
	CLOCKSTART;
	sort(reads->begin(),reads->end(), compareReads);	// Sort the reads lexicographically.
	CLOCKSTOP;
}

/**********************************************************************************************************************
	This function returns the number of reads
**********************************************************************************************************************/
UINT64 Dataset::getNumberOfReads(void)
{
	return numberOfReads;
}


/**********************************************************************************************************************
	This function returns the number of unique reads
**********************************************************************************************************************/
UINT64 Dataset::getNumberOfUniqueReads(void)
{
	return numberOfUniqueReads;
}

/**********************************************************************************************************************
	Returns true if the read contains only {A,C,G,T} and does not contain more than 80% of the same nucleotide
**********************************************************************************************************************/
bool Dataset::testRead(const string & read)
{

	UINT64 cnt[4] = {0,0,0,0};
	UINT64 readLength = read.length();
	for(UINT64 i = 0; i < readLength; i++) // Count the number of A's, C's , G's and T's in the string.
	{
		if(read[i]!= 'A' && read[i] != 'C' && read[i] != 'G' && read[i] != 'T')
			return false;
		cnt[(read[i] >> 1) & 0X03]++; // Least significant 2nd and 3rd bits of ASCII value used here
	}
	UINT64 threshold = read.length()*.8;	// 80% of the length.
	if(cnt[0] >= threshold || cnt[1] >= threshold || cnt[2] >= threshold || cnt[3] >= threshold)
		return false;	// If 80% bases are the same base.
	return true;
}




/**********************************************************************************************************************
	Search a read in the dataset using binary search
**********************************************************************************************************************/
Read * Dataset::getReadFromString(const string & read)
{
	UINT64 min = 0, max = getNumberOfUniqueReads()-1;
	string readReverse = reverseComplement(read);
	int comparator;
	if(read.compare(readReverse) < 0)
	{
		while (max >= min) 															// At first search for the forward string.
		{
			UINT64 mid = (min + max) / 2; 	// Determine which subarray to search.
			comparator = reads->at(mid)->getStringForward().compare(read.c_str());
			if(comparator == 0)
				return reads->at(mid);
			else if (comparator < 0) 	// Change min index to search upper subarray.
				min = mid + 1;
			else if (comparator > 0) 	// Change max index to search lower subarray.
				max = mid - 1;
		}
	}
	else
	{
		while (max >= min) 																	// If forward string is not found then search for the reverse string
		{
			UINT64 mid = (min+max) / 2; 													// Determine which subarray to search
			comparator = reads->at(mid)->getStringForward().compare(readReverse.c_str());
			if( comparator == 0)
				return reads->at(mid);
			else if (comparator < 0) 	// Change min index to search upper subarray.
				min = mid + 1;
			else if (comparator > 0) 	// Change max index to search lower subarray.
				max = mid - 1;
		}
	}
	MYEXIT("String not found in Dataset: "+read);
}




/**********************************************************************************************************************
	Returns the reverse complement of a read
**********************************************************************************************************************/
string Dataset::reverseComplement(const string & read)
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
	This function returns a read from its ID.
**********************************************************************************************************************/
Read * Dataset::getReadFromID(UINT64 ID)
{
	if(ID < 1 || ID > numberOfUniqueReads)	// ID outside the range.
	{
		stringstream ss;
		ss << "ID " << ID << " out of bound.";
		string s = ss.str();
		MYEXIT(s);
	}
	return reads->at(ID - 1);
}


/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Dataset::~Dataset(void)
{
	// Free the memory used by the dataset.
	for(UINT64 i = 0; i < reads->size(); i++)
	{
		delete reads->at(i);
	}
	delete reads;

}
