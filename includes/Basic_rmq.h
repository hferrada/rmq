#ifndef BASIC_DRF_H_
#define BASIC_DRF_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <dirent.h>

namespace rmqrmm {
	#ifndef uchar // 1 byte
	#define uchar unsigned char
	#endif
	#ifndef suint
	#define suint short unsigned int
	#endif
	#ifndef uint
	#define uint unsigned int
	#endif
	#ifndef ulong
	#define ulong unsigned long
	#endif

	const uint W64 = 64;
	const uint BW64 = 6;	// pow of two for W64
	const uint WW64 = 128;
	const uint W64minusone = 63;
	const uint W64m8 = 56;	// W64 - 8;
	const uint N8W64 = 8;	// W64/8 = 8
	const ulong maskW63 = 0x8000000000000000;

	// popcount array for uchars
	const unsigned char __popcount_tab[] = {
		0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
	};

	// reads i-th bit from e (left to right)
	#define readBit64(e,i) ((e[i/W64] >> (W64minusone-i%W64)) & 1)

	// sets bit i-ht in e (left to right)
	void setBit64(ulong *e, ulong i);

	// cleans bit i-ht in e (left to right)
	void cleanBit64(ulong *e, ulong i);

	// print the last cantBits of unisned int x
	void printBitsNum(uint x, uint cantBits);

	// print the last cantBits of unisned long int x
	void printBitsNum64(ulong x, uint cantBits);

	// print W64 bits of unsigned long int x
	void printBitsUlong(ulong x);

	// compute ceiling to logarithm base k for num
	uint ceilingLog64(ulong num, uint k);

	// return (in a unsigned integer) the number in A from bits of position 'ini' to 'ini+len-1'
	ulong getNum64(ulong *A, ulong ini, uint len);
	long int getNumLI64(long int *A, ulong ini, uint len);

	// Extract n cells: A[sp,...,sp+n-1] and stores values in B[0,...,n-1], where each cell has lenCell bits.
	void extractUlongs(ulong *A, ulong sp, ulong n, uint lenCell, ulong *B);

	// set the number x as a bitstring sequence in *A. In the range of bits [ini, .. ini+len-1] of *A. Here x has len bits
	void setNum64(ulong *A, ulong ini, uint len, ulong x);
	void setNumLI64(long int *A, ulong ini, uint len, long int x);

	// return in Milliseconds
	double getTime_ms();

	uint popcount_Rank32(uint x);
	uint popcount_Rank64(ulong x);

} /* namespace rmqrmm */

#endif /* BASIC_DRF_H_ */
