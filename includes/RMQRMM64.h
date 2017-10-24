/*
 * RMQRMM64.h
 *
 *  Created on: 16-06-2014
 *      Author: hector
 *
 * To create the BP sequence:
 * We use a stack Q to store positions of A during the construction, which will allocate
 * O(log n) cells in average. We first put a 0 value in Q for representing an extra cell A[0]=−∞
 * as the root of our virtual tree and put the first opening parenthesis in P (i.e., P [0] = 1).
 * Then for each cell A[k], k ≥ 1, we append a closing parenthesis for each value j>0 stored
 * in Q, such that A[j]≥A[k], deleting these values from Q. Next, we append an opening
 * parenthesis to represent A[k], store the k value in Q and continue with the next cell A[k+1]
 * if k<n. When we process all the cells of A, we append a closing parenthesis to P for each
 * value still stored in Q, and then we discard Q. As a result we have generated our balanced
 * parentheses sequence P[1..2n] in O(n) time, because we exactly insert/delete only once each
 * k value to/from Q, ∀ 1≤k≤ n.
*/

#ifndef RMQRMM64_H_
#define RMQRMM64_H_

#include "Basic_rmq.h"
using namespace std;
using namespace rmqrmm;

// ___________________________________________________
// Standard Configuration !
#define BLK 128		// size of blocks (BLK bits each one), (power of 2 >= W)
#define SS 128		// select_1 sampling size. Power of 2

// ___________________________________________________
// FIXIED VALUES:
#define BST 8		// bits for small block (for popcount)
#define BSTMOne 7	// BST minus one
// ___________________________________________________


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
// it stores the minimum in a samll block from right to left
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
	typedef struct StackBP{
		long int val;
		StackBP *next;
	} StackBP;

private:
	uint BBLK;				// bit of BLK = log(BLK)
	uint N8BLK;				// BLK/8;
	uint MBLK;				// BLK/2
	uint BMBLK;				// bit of BMBLK = log(BMBLK)
	uint BSS;				// bit of SS = log(SS)

	ulong nW;				// number of words to store P
	ulong *P;				// sequence of 2n' balanced parentheses, which represents a tree with n/2 nodes.

	ulong nBin;				// 0 < nBin <= n. This determine the partition of the interval for the Binary Tree and the last block
							// 'n - nBin' is the number of bits that must be considered as a single last block of length 'n - nBin'
	uint lenLast;			// length of lastBlock (number of bits)

	ulong cantN;			// total nodes = leaves + cantIN
	ulong cantIN;			// number of internal nodes of the min-max tree (nodes = cantIN + leaves)
	ulong leaves;			// number of leaves of min-max tree
	ulong leavesBottom;		// number of leaves in the last level h (perhaps there are leaves in the level h-1 too)
	ulong firstLeaf;		// position of the first leaf (the left-most leaf)

	ulong *BkM;				// the minimum excess value for internal nodes. MIN bits per element, where MIN is the logarithm of greater value
							// Bkwd_MinIN and Fwd_MinIN store the sign in the first bit for each item
	uint lgBkM;				// bit for each cell in BkM[]

	ulong nBLK;				// # BLK
	ulong *TMinB;			// Table of leaf minimum
	uint lg_MinB;
	ulong *TSumB;			// Table of global excess (semi sum) for each super block in P (groups of k blocks).
	uint lg_SumB;

	ulong *TSS;				// Sampling for select respect to Super Blocks. It stores number of SB instead of position of P.
	uint lenSS;				// # cells for TSS
	uint lg_SS;				// bit per cell in TSS = log(lenSS)

	uint sizeRMM;			// in bytes

public:
	ulong nP;				// Length of sequence P (n parentheses and n/2 nodes)

	static bool TRACE;		// true: print all details for console
	static bool RUNTEST;
	static bool SHOW_SIZE;
	static uint TEST;

	// Initializes variables
	void init(ulong len);

	// bitsPC are the number of bits for each cell in A[0..len-1]
	// if deleteA=true then the array A will the delete after to create the BP sequence P
	RMQRMM64(ulong *A, uint bitsPC, ulong len, bool deleteA);
	RMQRMM64(short int *A, ulong len);
	RMQRMM64(int *A, ulong len);
	RMQRMM64(long int *A, ulong len);
	RMQRMM64(char *fileName);
	virtual ~RMQRMM64();

	void createMinMaxTree();
	void createTables();

	// give the excess from 0 to pos
	long int sumAtPos(long int pos);
	void test_sumAtPos();

	ulong rank_1(ulong i);
	void test_rank_1();

	ulong select_1(ulong i);
	void test_select_1();

	// give the excess from block 0 to block blk
	long int sumAtBlock(long int blk);

	// return the position in the block 'blk' where is the minimum 'Min' of the block
	ulong positionMinblock(ulong blk);
	void test_positionMinblock();

	// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
	long int computeSumOfNode(ulong node, ulong dist);

	// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
	ulong leavesOfNode(ulong node, ulong dist, ulong *leafL);
	ulong computeLeavesOfNode(ulong node, ulong dist);

	// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
	void search_min_block(ulong x1, ulong x2, long int *min, long int *curSum, ulong *position);
	void test_search_min_block();

	// return the position of the open parenthesis closet to the root between i and j
	ulong rmqi(ulong i, ulong j);
	ulong rmqi_rmm(ulong x1, ulong x2, long int *min, long int *currSum, ulong posMin);
	void test_rmqi();

	// =================================================================================================================
	// query for RMQ
	ulong queryRMQ(ulong i, ulong j);

	uint getSize();

	// save the Data Structure in file 'fileName'
	void saveDS(char *fileName);

	// load the Data Structure from the file 'fileName'
	void loadDS(char *fileName);

	void printTree();
};

#endif /* RMQRMM64_H_ */
