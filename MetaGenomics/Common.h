/*
 * Common.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef COMMON_H_
#define COMMON_H_

#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <streambuf>
#include <sys/stat.h>
#include <map>
#include <string>
#include <cstring>
#include <hash_map>
#include <queue>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <functional>
#include <mpi.h>
#include <bitset>

using namespace std;

typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned long UINT32;
typedef long INT32;
typedef unsigned long long UINT64;
typedef long long INT64;


#define aStatisticsThreshold 3
#define minDelta 1000
#define deadEndLength 10
#define minimumSupport 3
#define loopLimit 15
#define	coverageDepthLB 40
#define coverageDepthUB 200
#define insertSizeRangeSD 3	// 3 means mean +/- 3 SD

//Multi-thread parallel options
#define DEF_THREAD_COUNT 4
#define MAX_PAR_GRAPH_SIZE 10000
#define MPI_BLOCK 10

//	Exit code that displays the place of exit and message.
#define MYEXIT(a) { cout << endl << "Exit from File: " << __FILE__ << " Line: " << __LINE__ << " Function: " << __FUNCTION__ << "()" << endl << "Message: " << a << endl; exit(0);}
// Print which function is currently executing. Only for functions that take long time


#define SSTR( x ) dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

// To keep time information of functions.
#define CLOCKSTART INT64 mem_start = checkMemoryUsage(); double begin = omp_get_wtime(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define CLOCKSTOP INT64 mem_end = checkMemoryUsage(); double end = omp_get_wtime(); cout << "Function " << __FUNCTION__ << "() finished in " <<(end - begin)<< " Seconds." << endl << "Memory used: " << mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl << endl;

// Get the memory usage with a Linux kernel.
inline unsigned int checkMemoryUsage()
{
    // get KB memory into count
    unsigned int count=0;

    #if defined(__linux__)
    ifstream f("/proc/self/status"); // read the linux file
    while(!f.eof()){
        string key;
        f>>key;
        if(key=="VmData:"){     // size of data
            f>>count;
        break;
        }
    }
    f.close();
    #endif

    // return MBs memory (size of data)
    return (count/1024);
};


inline void str_reverse( char *str ) {
    char *str_end = strchr( str, 0 );
    std::reverse( str, str_end );
}

inline std::string reverseComplement(const std::string & seq)
{
	uint64_t sLength = seq.length();
	std::string reverse(sLength,'0');
	for(uint64_t i = 0;i < sLength; i++)
	{
		if( seq[i] & 0X02 ) // C or G
			reverse.at(sLength -  i - 1 ) = seq[i] ^ 0X04;
		else // A <==> T
			reverse.at(sLength -  i - 1 ) = seq[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}

#endif /* COMMON_H_ */
