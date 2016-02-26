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
	superReadID = 0;
	read = new dna_bitset();
}

/**********************************************************************************************************************
	Another constructor
**********************************************************************************************************************/
Read::Read(const string & s)
{
	// Initialize the variables.
	read = new dna_bitset(s, s.length());
	readNumber = 0;
	readName="";
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
	dna_bitset *tmpRead = new dna_bitset(s, s.length());
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




