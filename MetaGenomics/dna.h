#ifndef __dna_h__
#define __dna_h__


#include <stdexcept>	/* std::invalid_argument */
#include <inttypes.h>
#include "Common.h"


#define BASE_MASK 0x00000003	/* binary: 11 */

#define FRAG_SIZE 32
#define HALF_FRAG_SIZE 16	/* Half of the number of bits used as the bit storage block. Here it is 8 so it is set to 4.*/

/* useful constants */
enum
{
	BASE_A = 0x0,	/* binary: 00 */
	BASE_C = 0x1,	/*'binary: 01 */
	BASE_G = 0x2,	/* binary: 10 */
	BASE_T = 0x3,	/* binary: 11 */
};

class dna_bitset
{
private:
	uint32_t* m_data;
	size_t   m_len;
public:

	/**
		 * @brief returns 1 is the substrings are equal
		 * @param this pointer to first sequence to be compared (subject)
		 * @param start index in subject to compare
		 * @param length of substring of subject to compare
		 * @param pointer to second sequence to be compared (query)
		 * @param start index in query to compare
		 * @param length of substring of query to compare
	*/
	bool compareSubString (UINT64 seq1Start, UINT64 seq1Len, const dna_bitset *seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient) const;
	/**
		 * @brief default constructor
	 */
	dna_bitset ()
	{
		m_len=0;
		m_data=NULL;
	}
	/**
	 * @brief constructor
	 * @param dna_str a string containing a DNA sequence (e.g. "ATGCA...")
	 * @param dna_len length of the DNA sequence
	 */
	dna_bitset (const string & dna_str, const size_t dna_len)
	{
		m_len = dna_len;

		/* number of bytes necessary to store dna_str as a bitset */
		size_t dna_word = (dna_len / HALF_FRAG_SIZE) + (dna_len % HALF_FRAG_SIZE != 0);

		m_data = new uint32_t[dna_word];

		std::memset(m_data, 0, dna_word*(FRAG_SIZE/8));

		/* for each base of the DNA sequence */
		for (size_t i = 0; i < dna_len; i++)
		{
			uint32_t shift = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);

			switch (dna_str[i])
			{
			case 'A':
				m_data[i/HALF_FRAG_SIZE] |= BASE_A << shift;
				break;
			case 'C':
				m_data[i/HALF_FRAG_SIZE] |= BASE_C << shift;
				break;
			case 'G':
				m_data[i/HALF_FRAG_SIZE] |= BASE_G << shift;
				break;
			case 'T':
				m_data[i/HALF_FRAG_SIZE] |= BASE_T << shift;
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}
		}
	}
	/**
	 * @brief destructor
	 */
	~dna_bitset ()
	{
		delete[] m_data;
	}
	/**
		 * @brief returns reverse complement of the stored DNA sequence as a C++ string
		 */
	string toRevComplement () const
	{
		string dna_str="";
		string strRevCArr[] = {"T","G","C","A"};
		/* for each base of the DNA sequence */
		for (int i = m_len-1; i >= 0; i--)
		{
			uint8_t shift = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);

			/* get the i-th DNA base */
			uint32_t base = (m_data[i/HALF_FRAG_SIZE] & (BASE_MASK << shift)) >> shift;
			dna_str += strRevCArr[base];
		}
		return dna_str;
	}
	/**
		* @brief returns the stored DNA sequence as a C++ string
	 */
	string toString () const
	{
		string dna_str="";
		string strArr[] = {"A","C","G","T"};
		/* for each base of the DNA sequence */
		for (size_t i = 0; i < m_len; i++)
		{
			uint8_t shift = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);
			/* get the i-th DNA base */
			uint32_t base = (m_data[i/HALF_FRAG_SIZE] & (BASE_MASK << shift)) >> shift;
			dna_str += strArr[base];
		}
		return dna_str;
	}
	int getLength(void) const {return m_len;}


};

#endif /* __dna_h__ */
