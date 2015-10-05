/*
 * RMQRMM64.h
 *
 *  Created on: 16-06-2014
 *      Author: hector
 */

#ifndef RMQRMM64_H_
#define RMQRMM64_H_

#include "Basic_rmq.h"
using namespace std;
using namespace rmqrmm;

#define Srmq 256	// size of blocks (s bits each one), (power of 2 >= W)
#define PotSrmq 8	// power for block = log(Srmq)
#define SrmqD 512	// 2*Srmq
#define SrmqM 128	// Srmq/2;
#define N8Srmq 32 	// Srmq/8;
#define SuBrmq 2	// number of leaves for each super block (power of 2 > RB)
#define BSrmq 8		// bits for each little block
#define BrmqMOne 7	// BSrmq minus one

const ulong RMMMasks[] = {0xFF00000000000000, 0x00FF000000000000, 0x0000FF0000000000, 0x000000FF00000000,
							 0x00000000FF000000, 0x0000000000FF0000, 0x000000000000FF00, 0x00000000000000FF,};

// -8 <= sum <= 8; // 2 bytes per cell --> 512 bytes
const short int T_SUM_BLOCK[] = {
		-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,
		-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
		-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
		-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8,
};

// -8 <= min <= 1; // 2 Bytes per cell --> 512 bytes
// it stores the relative level closest to the root for any 1, it no considers the 0's
const short int T_MIN_BCK[] = {
		0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,
		0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,
		0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,
		0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,
		0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,
		0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,
		0,1,1,3,1,3,3,5,1,3,3,5,3,5,5,7,
		0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,
		0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,
		0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,
		0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,
		0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
		0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8
};


// 8 bytes for each subset --> 256 x 8 = 2048 bytes
// bkwd_d[x][i] = position in x of (min(x)-i) = position in x of (minBkwd[x]-i), position from right to left, x has 8 bits
const uchar T_BCK_D[][8] = {
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,3,0,0,0,0,0,0,},{2,3,0,0,0,0,0,0,},{0,1,2,3,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{4,0,0,0,0,0,0,0,},{0,1,4,0,0,0,0,0,},{4,0,0,0,0,0,0,0,},{0,3,4,0,0,0,0,0,},{2,3,4,0,0,0,0,0,},{0,1,2,3,4,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,3,0,0,0,0,0,0,},{2,3,0,0,0,0,0,0,},{0,1,2,3,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,5,0,0,0,0,0,0,},{2,5,0,0,0,0,0,0,},{0,1,2,5,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,5,0,0,0,0,0,0,},{4,5,0,0,0,0,0,0,},{0,1,4,5,0,0,0,0,},{4,5,0,0,0,0,0,0,},{0,3,4,5,0,0,0,0,},{2,3,4,5,0,0,0,0,},{0,1,2,3,4,5,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,3,0,0,0,0,0,0,},{2,3,0,0,0,0,0,0,},{0,1,2,3,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{4,0,0,0,0,0,0,0,},{0,1,4,0,0,0,0,0,},{4,0,0,0,0,0,0,0,},{0,3,4,0,0,0,0,0,},{2,3,4,0,0,0,0,0,},{0,1,2,3,4,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{6,0,0,0,0,0,0,0,},{0,1,6,0,0,0,0,0,},{6,0,0,0,0,0,0,0,},{0,3,6,0,0,0,0,0,},{2,3,6,0,0,0,0,0,},{0,1,2,3,6,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{6,0,0,0,0,0,0,0,},{0,1,6,0,0,0,0,0,},{6,0,0,0,0,0,0,0,},{0,5,6,0,0,0,0,0,},{2,5,6,0,0,0,0,0,},{0,1,2,5,6,0,0,0,},
		{6,0,0,0,0,0,0,0,},{0,5,6,0,0,0,0,0,},{4,5,6,0,0,0,0,0,},{0,1,4,5,6,0,0,0,},{4,5,6,0,0,0,0,0,},{0,3,4,5,6,0,0,0,},{2,3,4,5,6,0,0,0,},{0,1,2,3,4,5,6,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,3,0,0,0,0,0,0,},{2,3,0,0,0,0,0,0,},{0,1,2,3,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{4,0,0,0,0,0,0,0,},{0,1,4,0,0,0,0,0,},{4,0,0,0,0,0,0,0,},{0,3,4,0,0,0,0,0,},{2,3,4,0,0,0,0,0,},{0,1,2,3,4,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,3,0,0,0,0,0,0,},{2,3,0,0,0,0,0,0,},{0,1,2,3,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,5,0,0,0,0,0,0,},{2,5,0,0,0,0,0,0,},{0,1,2,5,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,5,0,0,0,0,0,0,},{4,5,0,0,0,0,0,0,},{0,1,4,5,0,0,0,0,},{4,5,0,0,0,0,0,0,},{0,3,4,5,0,0,0,0,},{2,3,4,5,0,0,0,0,},{0,1,2,3,4,5,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{2,0,0,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,3,0,0,0,0,0,0,},{2,3,0,0,0,0,0,0,},{0,1,2,3,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,7,0,0,0,0,0,0,},{2,7,0,0,0,0,0,0,},{0,1,2,7,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,7,0,0,0,0,0,0,},{4,7,0,0,0,0,0,0,},{0,1,4,7,0,0,0,0,},{4,7,0,0,0,0,0,0,},{0,3,4,7,0,0,0,0,},{2,3,4,7,0,0,0,0,},{0,1,2,3,4,7,0,0,},
		{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,0,0,0,0,0,0,0,},{0,7,0,0,0,0,0,0,},{2,7,0,0,0,0,0,0,},{0,1,2,7,0,0,0,0,},
		{0,0,0,0,0,0,0,0,},{0,7,0,0,0,0,0,0,},{6,7,0,0,0,0,0,0,},{0,1,6,7,0,0,0,0,},{6,7,0,0,0,0,0,0,},{0,3,6,7,0,0,0,0,},{2,3,6,7,0,0,0,0,},{0,1,2,3,6,7,0,0,},
		{0,0,0,0,0,0,0,0,},{0,7,0,0,0,0,0,0,},{6,7,0,0,0,0,0,0,},{0,1,6,7,0,0,0,0,},{6,7,0,0,0,0,0,0,},{0,5,6,7,0,0,0,0,},{2,5,6,7,0,0,0,0,},{0,1,2,5,6,7,0,0,},
		{6,7,0,0,0,0,0,0,},{0,5,6,7,0,0,0,0,},{4,5,6,7,0,0,0,0,},{0,1,4,5,6,7,0,0,},{4,5,6,7,0,0,0,0,},{0,3,4,5,6,7,0,0,},{2,3,4,5,6,7,0,0,},{0,1,2,3,4,5,6,7,},
};

class RMQRMM64 {
private:
	ulong nW;				// number of words to store P
	ulong *P;				// sequence of 2n' balanced parentheses, which represents a tree with n/2 nodes.

	ulong rank1_Bin;		// number of 1's from 0 to nBin
	ulong nBin;				// 0 < nBin <= n. This determine the partition of the interval for the Binary Tree and the last block
							// 'n - nBin' is the number of bits that must be considered as a single last block of length 'n - nBin'

	uint lenLB;				// length of lastBlock (number of bits)
	uint bitsSuB, bitsRB;	// by default bitsSuB=256 and bitsRB=128
	uint h;					// nim-max tree's height

	ulong cantN;			// total nodes = leaves + cantIN
	ulong cantIN;			// number of internal nodes of the min-max tree (nodes = cantIN + leaves)
	ulong leaves;			// number of leaves of min-max tree
	ulong leavesBottom;		// number of leaves in the last level h (perhaps there are leaves in the level h-1 too)
	ulong firstLeaf;		// position of the first leaf (the left-most leaf)

	ulong *Bkwd_MinIN;		// the minimum excess value for internal nodes. MIN bits per element, where MIN is the logarithm of greater value
							// Bkwd_MinIN and Fwd_MinIN store the sign in the first bit for each item
	uint lgMIN_BCK;

	ulong *TSBlock;			// Table of global excess for each superblock in P (groups of k blocks). Size = (n/ks)*MAXEXC
	ulong lenSB;			// number of super/relative blocks
	uint MAX_SupB;
	uint lgMAX_SupB;

	char *TRBlock;			// Table of excess relative only for the first block in each super block.
							// TRBlock[i] = sum_i/2. where, sim_i is the relative sum for the first block of the superblock i,
							// also -255 <= sum_i <= 255 and sum_i is an even number. We need 8 bits for each value.
	ulong *Bfull;			// this indicates the leaves which contain in TRBlock the numbers -256 or 256
	uchar *TPMinB;			// Table of minimum positions for each leaf. values between 0 and 255
	ulong *TMinB;			// Table of leaf minimum

	uint MAX_B;				// the greater global excess for all block.
	uint lgMAX_B;
	uint MAX_BCK;			// the greater value for backward interval
	int MIN_BCK;			// the lowest value for backward interval

	uint sizeRMM;			// in bytes

public:
	ulong nP;				// Length of sequence P (n parentheses and n/2 nodes)

	static bool TRACE;		// true: print all details for console
	static bool RUNTEST;
	static uint TEST;

	RMQRMM64(long int *A, ulong len);
	RMQRMM64(char *fileName);
	virtual ~RMQRMM64();

	void createMinMaxTree();
	void createTables();
	void printTree();

	ulong binRank_1(ulong i);
	ulong binSelect_1(ulong i);

	ulong rank_1(ulong i);
	void test_rank_1();

	ulong select_1(ulong i);
	void test_select_1();

	// give the excess from 0 to pos
	long int sumAtPos(long int pos);
	void test_sumAtPos();

	// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
	long int computeSumOfNode(ulong node, ulong dist);

	// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
	ulong computeLeavesOfNode(ulong node, ulong dist);

	// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
	void search_min_block(ulong x1, ulong x2, long int *min, long int *curSum, ulong *position);

	// return the position of the open parenthesis closet to the root between i and j
	ulong rmqi(ulong i, ulong j);
	ulong rmqi_rmm(ulong x1, ulong x2, long int *min, long int *currSum, ulong posMin);
	void test_rmqi();

	// query for RMQ
	ulong queryRMQ(ulong i, ulong j);

	uint getSize();

	// save the Data Structure in file 'fileName'
	void saveDS(char *fileName);

	// load the Data Structure from the file 'fileName'
	void loadDS(char *fileName);
};

#endif /* RMQRMM64_H_ */
