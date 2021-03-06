/*
 * RangeMMTree64.h
 *
 *  Created on: 01-07-2014
 *      Author: hector
 */

#ifndef RANGEMMTREE64_H_
#define RANGEMMTREE64_H_

#include "Basic_drf64.h"

using namespace std;
using namespace drf64;

#define BitB 8		// bits for each little block
#define BmOne 7		// BitB minus one
#define S 128		// size of blocks (s bits each one)
#define PotS 7		// power for block = log(S)
#define N8S 16		// N8S = S/8 = 4

const ulong RMM_Masks[] = {0xFF00000000000000, 0x00FF000000000000, 0x0000FF0000000000, 0x000000FF00000000,
						0x00000000FF000000, 0x0000000000FF0000, 0x000000000000FF00, 0x00000000000000FF,};

// length=256; -8 <= min <= 1; // 2 bytes per cell --> 512 bytes
const uchar T_MIN_FWDI[] = {
		8,7,6,6,6,5,5,5,6,5,4,4,4,4,4,4,
		6,5,4,4,4,3,3,3,4,3,3,3,3,3,3,3,
		6,5,4,4,4,3,3,3,4,3,2,2,2,2,2,2,
		4,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
		6,5,4,4,4,3,3,3,4,3,2,2,2,2,2,2,
		4,3,2,2,2,1,1,1,2,1,1,1,1,1,1,1,
		4,3,2,2,2,1,1,1,2,1,1,1,1,1,1,1,
		2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
		6,5,4,4,4,3,3,3,4,3,2,2,2,2,2,2,
		4,3,2,2,2,1,1,1,2,1,1,1,1,1,1,1,
		4,3,2,2,2,1,1,1,2,1,0,0,0,0,0,0,
		2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		4,3,2,2,2,1,1,1,2,1,0,0,0,0,0,0,
		2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

// length=256; -8 <= sum <= 8; // 2 bytes per cell --> 512 bytes
const short int T_SUM_BLOCKI[] = {
		8,6,6,4,6,4,4,2,6,4,4,2,4,2,2,0,
		6,4,4,2,4,2,2,0,4,2,2,0,2,0,0,-2,
		6,4,4,2,4,2,2,0,4,2,2,0,2,0,0,-2,
		4,2,2,0,2,0,0,-2,2,0,0,-2,0,-2,-2,-4,
		6,4,4,2,4,2,2,0,4,2,2,0,2,0,0,-2,
		4,2,2,0,2,0,0,-2,2,0,0,-2,0,-2,-2,-4,
		4,2,2,0,2,0,0,-2,2,0,0,-2,0,-2,-2,-4,
		2,0,0,-2,0,-2,-2,-4,0,-2,-2,-4,-2,-4,-4,-6,
		6,4,4,2,4,2,2,0,4,2,2,0,2,0,0,-2,
		4,2,2,0,2,0,0,-2,2,0,0,-2,0,-2,-2,-4,
		4,2,2,0,2,0,0,-2,2,0,0,-2,0,-2,-2,-4,
		2,0,0,-2,0,-2,-2,-4,0,-2,-2,-4,-2,-4,-4,-6,
		4,2,2,0,2,0,0,-2,2,0,0,-2,0,-2,-2,-4,
		2,0,0,-2,0,-2,-2,-4,0,-2,-2,-4,-2,-4,-4,-6,
		2,0,0,-2,0,-2,-2,-4,0,-2,-2,-4,-2,-4,-4,-6,
		0,-2,-2,-4,-2,-4,-4,-6,-2,-4,-4,-6,-4,-6,-6,-8,
};

// 8 bytes for each subset --> 256 x 8 = 2048 bytes
// length=8*256; fwd_d[x][i] = position in x of (min(x)+i) = position in x of (minFwd[x]+i), position from left to right, x has 8 bits
const uchar T_FWD_D[][8] = {
		{7,6,5,4,3,2,1,0,},{6,5,4,3,2,1,0,0,},{5,4,3,2,1,0,0,0,},{5,4,3,2,1,0,0,0,},{7,4,3,2,1,0,0,0,},{4,3,2,1,0,0,0,0,},{4,3,2,1,0,0,0,0,},{4,3,2,1,0,0,0,0,},
		{7,6,3,2,1,0,0,0,},{6,3,2,1,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,7,0,0,0,},
		{7,6,5,2,1,0,0,0,},{6,5,2,1,0,0,0,0,},{5,2,1,0,0,0,0,0,},{5,2,1,0,0,0,0,0,},{7,2,1,0,0,0,0,0,},{2,1,0,0,0,0,0,0,},{2,1,0,0,0,0,0,0,},{2,1,0,7,0,0,0,0,},
		{7,2,1,0,0,0,0,0,},{2,1,0,0,0,0,0,0,},{2,1,0,0,0,0,0,0,},{2,1,0,7,0,0,0,0,},{2,1,0,5,0,0,0,0,},{2,1,0,5,0,0,0,0,},{2,1,0,5,6,0,0,0,},{2,1,0,5,6,7,0,0,},
		{7,6,5,4,1,0,0,0,},{6,5,4,1,0,0,0,0,},{5,4,1,0,0,0,0,0,},{5,4,1,0,0,0,0,0,},{7,4,1,0,0,0,0,0,},{4,1,0,0,0,0,0,0,},{4,1,0,0,0,0,0,0,},{4,1,0,7,0,0,0,0,},
		{7,6,1,0,0,0,0,0,},{6,1,0,0,0,0,0,0,},{1,0,0,0,0,0,0,0,},{1,0,7,0,0,0,0,0,},{1,0,5,0,0,0,0,0,},{1,0,5,0,0,0,0,0,},{1,0,5,6,0,0,0,0,},{1,0,5,6,7,0,0,0,},
		{7,6,1,0,3,0,0,0,},{6,1,0,3,0,0,0,0,},{1,0,3,0,0,0,0,0,},{1,0,3,0,0,0,0,0,},{1,0,3,0,0,0,0,0,},{1,0,3,0,0,0,0,0,},{1,0,3,6,0,0,0,0,},{1,0,3,6,7,0,0,0,},
		{1,0,3,4,0,0,0,0,},{1,0,3,4,0,0,0,0,},{1,0,3,4,0,0,0,0,},{1,0,3,4,7,0,0,0,},{1,0,3,4,5,0,0,0,},{1,0,3,4,5,0,0,0,},{1,0,3,4,5,6,0,0,},{1,0,3,4,5,6,7,0,},
		{7,6,5,4,3,0,1,0,},{6,5,4,3,0,1,0,0,},{5,4,3,0,1,0,0,0,},{5,4,3,0,1,0,0,0,},{7,4,3,0,1,0,0,0,},{4,3,0,1,0,0,0,0,},{4,3,0,1,0,0,0,0,},{4,3,0,1,0,0,0,0,},
		{7,6,3,0,1,0,0,0,},{6,3,0,1,0,0,0,0,},{3,0,1,0,0,0,0,0,},{3,0,1,0,0,0,0,0,},{3,0,1,0,0,0,0,0,},{3,0,1,0,0,0,0,0,},{3,0,1,6,0,0,0,0,},{3,0,1,6,7,0,0,0,},
		{7,6,5,0,1,0,0,0,},{6,5,0,1,0,0,0,0,},{5,0,1,0,0,0,0,0,},{5,0,1,0,0,0,0,0,},{7,0,1,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,1,6,0,0,0,0,0,},{0,1,6,7,0,0,0,0,},
		{7,0,1,4,0,0,0,0,},{0,1,4,0,0,0,0,0,},{0,1,4,0,0,0,0,0,},{0,1,4,7,0,0,0,0,},{0,1,4,5,0,0,0,0,},{0,1,4,5,0,0,0,0,},{0,1,4,5,6,0,0,0,},{0,1,4,5,6,7,0,0,},
		{7,6,5,0,1,2,0,0,},{6,5,0,1,2,0,0,0,},{5,0,1,2,0,0,0,0,},{5,0,1,2,0,0,0,0,},{7,0,1,2,0,0,0,0,},{0,1,2,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},{0,1,2,7,0,0,0,0,},
		{7,0,1,2,0,0,0,0,},{0,1,2,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},{0,1,2,7,0,0,0,0,},{0,1,2,5,0,0,0,0,},{0,1,2,5,0,0,0,0,},{0,1,2,5,6,0,0,0,},{0,1,2,5,6,7,0,0,},
		{7,0,1,2,3,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,6,0,0,0,},{0,1,2,3,6,7,0,0,},
		{0,1,2,3,4,0,0,0,},{0,1,2,3,4,0,0,0,},{0,1,2,3,4,0,0,0,},{0,1,2,3,4,7,0,0,},{0,1,2,3,4,5,0,0,},{0,1,2,3,4,5,0,0,},{0,1,2,3,4,5,6,0,},{0,1,2,3,4,5,6,7,},
		{7,6,5,4,3,2,1,0,},{6,5,4,3,2,1,0,0,},{5,4,3,2,1,0,0,0,},{5,4,3,2,1,0,0,0,},{7,4,3,2,1,0,0,0,},{4,3,2,1,0,0,0,0,},{4,3,2,1,0,0,0,0,},{4,3,2,1,0,0,0,0,},
		{7,6,3,2,1,0,0,0,},{6,3,2,1,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,0,0,0,0,},{3,2,1,0,7,0,0,0,},
		{7,6,5,2,1,0,0,0,},{6,5,2,1,0,0,0,0,},{5,2,1,0,0,0,0,0,},{5,2,1,0,0,0,0,0,},{7,2,1,0,0,0,0,0,},{2,1,0,0,0,0,0,0,},{2,1,0,0,0,0,0,0,},{2,1,0,7,0,0,0,0,},
		{7,2,1,0,0,0,0,0,},{2,1,0,0,0,0,0,0,},{2,1,0,0,0,0,0,0,},{2,1,0,7,0,0,0,0,},{2,1,0,5,0,0,0,0,},{2,1,0,5,0,0,0,0,},{2,1,0,5,6,0,0,0,},{2,1,0,5,6,7,0,0,},
		{7,6,5,4,1,0,0,0,},{6,5,4,1,0,0,0,0,},{5,4,1,0,0,0,0,0,},{5,4,1,0,0,0,0,0,},{7,4,1,0,0,0,0,0,},{4,1,0,0,0,0,0,0,},{4,1,0,0,0,0,0,0,},{4,1,0,7,0,0,0,0,},
		{7,6,1,0,0,0,0,0,},{6,1,0,0,0,0,0,0,},{1,0,0,0,0,0,0,0,},{1,0,7,0,0,0,0,0,},{1,0,5,0,0,0,0,0,},{1,0,5,0,0,0,0,0,},{1,0,5,6,0,0,0,0,},{1,0,5,6,7,0,0,0,},
		{7,6,1,0,3,0,0,0,},{6,1,0,3,0,0,0,0,},{1,0,3,0,0,0,0,0,},{1,0,3,0,0,0,0,0,},{1,0,3,0,0,0,0,0,},{1,0,3,0,0,0,0,0,},{1,0,3,6,0,0,0,0,},{1,0,3,6,7,0,0,0,},
		{1,0,3,4,0,0,0,0,},{1,0,3,4,0,0,0,0,},{1,0,3,4,0,0,0,0,},{1,0,3,4,7,0,0,0,},{1,0,3,4,5,0,0,0,},{1,0,3,4,5,0,0,0,},{1,0,3,4,5,6,0,0,},{1,0,3,4,5,6,7,0,},
		{7,6,5,4,3,0,1,0,},{6,5,4,3,0,1,0,0,},{5,4,3,0,1,0,0,0,},{5,4,3,0,1,0,0,0,},{7,4,3,0,1,0,0,0,},{4,3,0,1,0,0,0,0,},{4,3,0,1,0,0,0,0,},{4,3,0,1,0,0,0,0,},
		{7,6,3,0,1,0,0,0,},{6,3,0,1,0,0,0,0,},{3,0,1,0,0,0,0,0,},{3,0,1,0,0,0,0,0,},{3,0,1,0,0,0,0,0,},{3,0,1,0,0,0,0,0,},{3,0,1,6,0,0,0,0,},{3,0,1,6,7,0,0,0,},
		{7,6,5,0,1,0,0,0,},{6,5,0,1,0,0,0,0,},{5,0,1,0,0,0,0,0,},{5,0,1,0,0,0,0,0,},{7,0,1,0,0,0,0,0,},{0,1,0,0,0,0,0,0,},{0,1,6,0,0,0,0,0,},{0,1,6,7,0,0,0,0,},
		{7,0,1,4,0,0,0,0,},{0,1,4,0,0,0,0,0,},{0,1,4,0,0,0,0,0,},{0,1,4,7,0,0,0,0,},{0,1,4,5,0,0,0,0,},{0,1,4,5,0,0,0,0,},{0,1,4,5,6,0,0,0,},{0,1,4,5,6,7,0,0,},
		{7,6,5,0,1,2,0,0,},{6,5,0,1,2,0,0,0,},{5,0,1,2,0,0,0,0,},{5,0,1,2,0,0,0,0,},{7,0,1,2,0,0,0,0,},{0,1,2,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},{0,1,2,7,0,0,0,0,},
		{7,0,1,2,0,0,0,0,},{0,1,2,0,0,0,0,0,},{0,1,2,0,0,0,0,0,},{0,1,2,7,0,0,0,0,},{0,1,2,5,0,0,0,0,},{0,1,2,5,0,0,0,0,},{0,1,2,5,6,0,0,0,},{0,1,2,5,6,7,0,0,},
		{7,0,1,2,3,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,0,0,0,0,},{0,1,2,3,6,0,0,0,},{0,1,2,3,6,7,0,0,},
		{0,1,2,3,4,0,0,0,},{0,1,2,3,4,0,0,0,},{0,1,2,3,4,0,0,0,},{0,1,2,3,4,7,0,0,},{0,1,2,3,4,5,0,0,},{0,1,2,3,4,5,0,0,},{0,1,2,3,4,5,6,0,},{0,1,2,3,4,5,6,7,},
};

// prev array for uchars. Size 256*8 bits
const uchar prev_tab[] = {
	0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
	6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
	8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
};


class RangeMMTree64 {
public:
	ulong *P;				// The length of the DFUDS sequences P is n, which represents a tree with n/2 nodes.
	uchar *labels;			// labels sequence
	ulong n;				// Length of sequence P (n parentheses and n/2 nodes)
	ulong nW;				// number of words to store P

	ulong rank1_Bin;		// number of 1's from 0 to nBin
	ulong nBin;				// 0 < nBin <= n. This determine the partition of the interval for the Binary Tree and the last block
							// 'n - nBin' is the number of bits that must be considered as a single last block of length 'n - nBin'
	uint lenLB;			// length of lastBlock (number of bits)

	// for nBin... first bits
	uint h;					// nim-max tree's height
	ulong cantN;			// total nodes = leaves + cantIN
	ulong cantIN;			// number of internal nodes of the min-max tree (nodes = cantIN + leaves)
	ulong leaves;			// number of leaves of min-max tree
	ulong leavesBottom;		// number of leaves in the last level h (perhaps there are leaves in the level h-1 too)
	ulong firstLeaf;		// position of the first leaf (the left-most leaf)

	ulong *Fwd_MinIN;		// the minimum excess value for internal nodes. MIN bits per element, where MIN is the logarithm of greater value
	uint lgMIN_FWD;

	ulong *TRBlock;			// Table of excess relative of each super block in P. Size = (n/s) log(2s+1)
	ulong *TSBlock;			// Table of global excess for each superblock in P (groups of k blocks). Size = (n/ks)*MAXEXC
	ulong lenSB;		// number of super/relative blocks
	ulong MAX_RelB;			// the greater relative value for TRBlock.
	uint lgMAX_RelB;
	ulong MAX_SupB;
	uint lgMAX_SupB;

	ulong sizeRMM;			// in bytes

	static bool TRACE;		// true: print all details for console
	static bool RUNTEST;
	static uint TEST;

	RangeMMTree64(ulong *PSeq, ulong len, uchar *labelSeq, bool showSize);
	RangeMMTree64(char *fileTrie, bool showSize);
	void createMinMaxTree(bool showSize);
	void createTables(bool showSize);

	// give the excess from 0 to pos
	long int sumAtPos(ulong pos);
	void test_sumAtPos();

	ulong binRank_1(ulong i);
	ulong binRank_0(ulong i);
	ulong rank_1(ulong i);
	void test_rank_1();

	ulong binSelect_0(ulong i);
	ulong select_0(ulong i);
	void test_select_0();

	// return the position >= k that is a close parenthesis (or 0)
	ulong selectNext0(const ulong k) const;
	void test_selectNext0();

	bool binFwd_search(ulong ini, ulong *d, ulong* pos);
	bool fwd_search(ulong x, ulong *d, ulong *pos);
	void test_fwd_search();


	// It gives the sum for this node, 'node=preorder+1'
	long int sumOfNode(ulong node);

	// return true and the position 'pos', pos < x+len, where the sum from x to x+len is d. If not found then return false
	bool fwd_block(ulong x, ulong *d, ulong* pos);
	bool fwd_Lastblock(ulong x, ulong d, ulong *pos);

	// return true and the relative position of the child (in 'child') labeled with 'c' for the node 'x'. The parameter 'end' is the position of the last child of x
	bool thereIsChild(ulong x, uchar c, ulong *child, uint len);

	// return true if this fictitious node has a unique child
	bool hasUniqueChild(ulong x);

	// it gives the size of the subtree rooted at x (without fictitious nodes)
	ulong subTreeSize(ulong x);

	void printTree();

	// save the Data Structure in file 'fileName' and return the amount of bytes saved
	ulong saveDS(char *fileName, bool showSize);

	// load the Data Structure from the file 'fileName'
	void loadDS(char *fileName, bool showSize);

	virtual ~RangeMMTree64();
};

#endif /* RANGEMMTREE64_H_ */
