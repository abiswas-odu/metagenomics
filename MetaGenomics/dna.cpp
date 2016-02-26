

#include "dna.h"

/**
	 * @brief returns 1 is the substrings are equal
	 * @param pointer to first sequence to be compared (subject)
	 * @param start index in subject to compare
	 * @param length of substring of subject to compare
	 * @param pointer to second sequence to be compared (query)
	 * @param start index in query to compare
	 * @param length of substring of query to compare
	 * @param overlap orientation
	 * // orient 0
	//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
	//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
	//				OR
	// orient 2
	//	 >---*****MMMMMMMMMMMMMMM*******------> read1
	//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of read2
	 *
	 * // orient 1
	//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
	//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
	//				OR
	// orient 3
	//	 >---*****MMMMMMMMMMMMMMM-------------> read1
	//		<*****MMMMMMMMMMMMMMM				Reverse Complement of Read2
*/
/*bool dna_bitset::compareSubString (UINT64 seq1Start, UINT64 seq1Len, const dna_bitset *seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient) const
{
	if(seq1Len!=seq2Len)
		return false;
	if(orient == 0 || orient== 1)
	{
		for (int i = seq1Start, j = seq2Start; i < seq1Start+seq1Len; i++,j++)
		{
			uint8_t shift1 = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);
			///* get the i-th DNA base from subject
			uint8_t base1 = (m_data[i/HALF_FRAG_SIZE] & (BASE_MASK << shift1)) >> shift1;

			uint8_t shift2 = (HALF_FRAG_SIZE*2-2) - 2*(j % HALF_FRAG_SIZE);
			///* get the j-th DNA base from query
			uint8_t base2 = (seq2->m_data[j/HALF_FRAG_SIZE] & (BASE_MASK << shift2)) >> shift2;

			if(base1!=base2)
				return false;
		}
	}
	else
	{
		for (int i = seq1Start, j = seq2->m_len-seq2Start-1; i < seq1Start+seq1Len; i++,j--)
		{
			uint8_t shift1 = (HALF_FRAG_SIZE*2-2) - 2*(i % HALF_FRAG_SIZE);
			///* get the i-th DNA base from subject
			uint8_t base1 = (m_data[i/HALF_FRAG_SIZE] & (BASE_MASK << shift1)) >> shift1;

			uint8_t shift2 = (HALF_FRAG_SIZE*2-2) - 2*(j % HALF_FRAG_SIZE);
			///* get the j-th DNA base from query
			uint8_t base2 = ~((seq2->m_data[j/HALF_FRAG_SIZE] & (BASE_MASK << shift2)) >> shift2) & BASE_MASK;

			if(base1!=base2)
				return false;
		}
	}
	return true;
}*/

bool dna_bitset::compareSubString (UINT64 seq1Start, UINT64 seq1Len, const dna_bitset *seq2, UINT64 seq2Start, UINT64 seq2Len, UINT64 orient) const
{

	if(seq1Len!=seq2Len)
		return false;
	if(orient == 0 || orient== 1)
	{
		UINT64 seq1Ei = seq1Start+seq1Len-1;
		UINT64 seq2Ei = seq2Start+seq2Len-1;
		for (size_t i = seq1Start, j = seq2Start; i <= seq1Ei; i+=HALF_FRAG_SIZE,j+=HALF_FRAG_SIZE)
		{
			uint8_t baseShift = i%HALF_FRAG_SIZE;
			uint32_t X1 = (m_data[i/HALF_FRAG_SIZE]) << 2*baseShift;				//shift bits to beginning of word boundary   00001111 --> 11110000
			uint32_t X2 = 0x00000000;
			if((i+(HALF_FRAG_SIZE-baseShift)) <= seq1Ei)							//move bits from next word into this word based on the number needed
				X2 = (m_data[i/HALF_FRAG_SIZE +1]) >> (FRAG_SIZE - 2*baseShift);

			uint32_t comp1 = X1 | X2;												//OR the two words to get the DNA bits of word length

			int maskEIndx = (i+HALF_FRAG_SIZE-1) - seq1Ei;				//If end index was reached, mask extra bits from that index onwards
			if(maskEIndx > 0)
				comp1 = comp1 & (0xFFFFFF << (uint32_t)maskEIndx*2);

			baseShift = j%HALF_FRAG_SIZE;
			X1 = (seq2->m_data[j/HALF_FRAG_SIZE]) << 2*baseShift;
			X2 = 0x00000000;
			if((j+(HALF_FRAG_SIZE-baseShift)) <= seq2Ei)
				X2 = (seq2->m_data[j/HALF_FRAG_SIZE +1]) >> (FRAG_SIZE - 2*baseShift);

			uint32_t comp2 = X1 | X2;

			maskEIndx = (j+HALF_FRAG_SIZE-1) - seq2Ei;
			if(maskEIndx > 0)
				comp2 = comp2 & (0xFFFFFFFF << (uint32_t)maskEIndx*2);

			if(comp1^comp2)
				return false;
		}
	}
	else
	{
		UINT64 seq1Ei = seq1Start+seq1Len-1;
		UINT64 seq2Ei = seq2->m_len-(seq2Start+seq2Len);
		for (size_t i = seq1Start, j = seq2->m_len-seq2Start-1; i <= seq1Ei; i+=HALF_FRAG_SIZE,j-=HALF_FRAG_SIZE)
		{
			uint8_t baseShift = i%HALF_FRAG_SIZE;
			uint32_t X1 = (m_data[i/HALF_FRAG_SIZE]) << 2*baseShift;
			uint32_t X2 = 0x00000000;
			if((i+(HALF_FRAG_SIZE-baseShift)) <= seq1Ei)
				X2 = (m_data[i/HALF_FRAG_SIZE +1]) >> (FRAG_SIZE - 2*baseShift);

			uint32_t comp1 = X1 | X2;

			int maskEIndx = (i+HALF_FRAG_SIZE-1) - seq1Ei;				//If end index was reached, mask extra bits from that index onwards
			if(maskEIndx > 0)
				comp1 = comp1 & (0xFFFFFF << (uint32_t)maskEIndx*2);

			baseShift = (j%HALF_FRAG_SIZE +1);
			X1 = (seq2->m_data[j/HALF_FRAG_SIZE]) >> (FRAG_SIZE - baseShift*2);
			X2 = 0x00000000;
			if((j-baseShift) >= seq2Ei)
				X2 = (seq2->m_data[j/HALF_FRAG_SIZE - 1]) << (baseShift*2);

			uint32_t comp2 = ~(X1 | X2);

			maskEIndx = seq2Ei - (j-HALF_FRAG_SIZE + 1);
			if(maskEIndx > 0)
				comp2 = comp2 & (0xFFFFFFFF >> (uint32_t)maskEIndx*2);
			//Reverse 32 bit value
			//comp2 = ((comp2 & 0x55555555) << 1) | ((comp2 & 0xAAAAAAAA) >> 1);
			comp2 = ((comp2 & 0x33333333) << 2) | ((comp2 & 0xCCCCCCCC) >> 2);
			comp2 = ((comp2 & 0x0F0F0F0F) << 4) | ((comp2 & 0xF0F0F0F0) >> 4);
			comp2 = ((comp2 & 0x00FF00FF) << 8) | ((comp2 & 0xFF00FF00) >> 8);
			comp2 = ((comp2 & 0x0000FFFF) << 16) | ((comp2 & 0xFFFF0000) >> 16);

			if(comp1^comp2)
				return false;
		}
	}
	return true;
}
