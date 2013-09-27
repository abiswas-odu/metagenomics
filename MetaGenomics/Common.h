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

using namespace std;

typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned long UINT32;
typedef long INT32;
typedef unsigned long long UINT64;
typedef long long INT64;


#define aStatisticsThreshold 3  // This in not used. Used previously in estimating the genome size in isolated genome assembly.
#define minDelta 1000			// This in not used. Used previously in estimating the genome size in isolated genome assembly.
#define deadEndLength 10		// Threshold of number of nodes in an edge for deleting dead-ends.
								// CP: if a node has many duplicate reads, it's still counted as one, right?
								// BH: Yes we only consider a read once.
#define minimumSupport 1		// Minimum number of support required to merge to edges (both using matepair paths and scaffolding step)
								// CP: should we consider using different thresholds for matepair paths and scaffolding?
								// BH: If we want we can use different thresholds.
#define loopLimit 15			// Number of time to loop in the main function.
#define insertSizeRangeSD 2		// 2 means mean +/- 2 SD
#define	coverageDepthLB 40		// Lower bound of coverage depth for the genome that we want to assemble.
#define coverageDepthUB 200	// Upper bound of the coverage depth of teh genome that we want to assemble.

//	Exit code that displays the place of exit and message.
#define MYEXIT(a) { cout << endl << "Exit from File: " << __FILE__ << " Line: " << __LINE__ << " Function: " << __FUNCTION__ << "()" << endl << "Message: " << a << endl; exit(0);}
// Print which function is currently executing. Only for functions that take long time


// To keep time information of functions.
#define CLOCKSTART INT64 mem_start = checkMemoryUsage(); clock_t begin = clock(); cout<<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()" << endl;
#define CLOCKSTOP INT64 mem_end = checkMemoryUsage(); clock_t end = clock(); cout << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) / CLOCKS_PER_SEC<< " Seconds." << endl << "Memory used: " << mem_end << " - " <<  mem_start << " = "<< mem_end - mem_start << " MB."<< endl << endl;

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



#endif /* COMMON_H_ */
