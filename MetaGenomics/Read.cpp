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
	isContainedRead = false;
	superReadID = 0;
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
	isContainedRead = false;
	superReadID = 0;
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Read::~Read(void)
{
	// delete all the pointers.
	delete read;
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
bool Read::compareReadOverlap(UINT64 seq1Start, UINT64 seq1Len, Read * seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient)
{
	return read->compareSubString(seq1Start,seq1Len,seq2->read,seq2Start,seq2Len, orient);
}




