#ifndef __dna_h__
#define __dna_h__


#include <stdexcept>	/* std::invalid_argument */
#include <inttypes.h>
#include "Common.h"


#define BASE_MASK 0x3	/* binary: 11 */

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
	uint8_t* m_data;
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
	dna_bitset (const char* dna_str, const size_t dna_len)
	{
		m_len = dna_len;

		/* number of bytes necessary to store dna_str as a bitset */
		size_t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);

		m_data = new uint8_t[dna_bytes];

		std::memset(m_data, 0, dna_bytes);

		/* for each base of the DNA sequence */
		for (size_t i = 0; i < dna_len; i++)
		{
			uint8_t shift = 6 - 2*(i % 4);

			switch (dna_str[i])
			{
			case 'A':
				m_data[i/4] |= BASE_A << shift;
				break;
			case 'C':
				m_data[i/4] |= BASE_C << shift;
				break;
			case 'G':
				m_data[i/4] |= BASE_G << shift;
				break;
			case 'T':
				m_data[i/4] |= BASE_T << shift;
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}

			shift = (shift == 0) ? 6 : shift-2;
		}
	}


	/**
	 * @brief setter method
	 * @param dna_str a string containing a DNA sequence (e.g. "ATGCA...")
	 * @param dna_len length of the DNA sequence
	 */
	void resetRead(const char* dna_str, const size_t dna_len)
	{
		m_len = dna_len;

		/* number of bytes necessary to store dna_str as a bitset */
		size_t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);

		m_data = new uint8_t[dna_bytes];

		std::memset(m_data, 0, dna_bytes);

		/* for each base of the DNA sequence */
		for (size_t i = 0; i < dna_len; i++)
		{
			uint8_t shift = 6 - 2*(i % 4);

			switch (dna_str[i])
			{
			case 'A':
				m_data[i/4] |= BASE_A << shift;
				break;
			case 'C':
				m_data[i/4] |= BASE_C << shift;
				break;
			case 'G':
				m_data[i/4] |= BASE_G << shift;
				break;
			case 'T':
				m_data[i/4] |= BASE_T << shift;
				break;
			default:
				throw std::invalid_argument("invalid DNA base");
			}

			shift = (shift == 0) ? 6 : shift-2;
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
			uint8_t shift = 6 - 2*(i % 4);

			/* get the i-th DNA base */
			uint8_t base = (m_data[i/4] & (BASE_MASK << shift)) >> shift;
			dna_str += strRevCArr[base];
		}
		return dna_str;
	}

	/**
	 * @brief returns the stored DNA sequence as a C string
	 */
	char* to_c_string () const
	{
		char* dna_str = new char[m_len+1];
		char *ch;
		/* for each base of the DNA sequence */
		for (size_t i = 0; i < m_len; i++)
		{
			uint8_t shift = 6 - 2*(i % 4);

			/* get the i-th DNA base */
			uint8_t base = (m_data[i/4] & (BASE_MASK << shift)) >> shift;

			switch (base)
			{
			case BASE_A:
				dna_str[i] = 'A';
				break;
			case BASE_C:
				dna_str[i] = 'C';
				break;
			case BASE_G:
				dna_str[i] = 'G';
				break;
			case BASE_T:
				dna_str[i] = 'T';
				break;
			default:
				throw std::runtime_error("Invalid DNA base");
			}
		}
		dna_str[m_len] = '\0';
		str_reverse(dna_str);
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
			uint8_t shift = 6 - 2*(i % 4);
			/* get the i-th DNA base */
			uint8_t base = (m_data[i/4] & (BASE_MASK << shift)) >> shift;
			dna_str += strArr[base];
		}
		return dna_str;
	}
	int getLength(void) const {return m_len;}


};

#endif /* __dna_h__ */
