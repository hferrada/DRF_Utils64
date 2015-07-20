/*
 * RangeMMTree64.cpp
 *
 *  Created on: 01-07-2014
 *      Author: hector
 */

#include "includes/RangeMMTree64.h"
bool RangeMMTree64::TRACE = false;
bool RangeMMTree64::RUNTEST = false;
uint RangeMMTree64::TEST = 1000;

RangeMMTree64::RangeMMTree64(char *fileTrie, bool showSize){
	loadDS(fileTrie, showSize);

	if (TRACE){
		ulong nb = n / S;
		if (nb%2 && nb > 1)
			nb--;
		printTree();
		cout << "# blocks " << nb << endl;
		cout << "rank1_Bin (1's) " << rank1_Bin << endl;
		cout << "last bit of binary tree: " << nBin << endl;
		cout << "Length of the last block " << lenLB << endl;
		cout << "Sum relative to binary tree: " << sumAtPos(nBin) << endl;
		cout << "======================"<< endl;
	}

	if (RUNTEST){
		cout << "running test... " << endl;
		test_sumAtPos();
		test_select_0();
		test_selectNext0();
		test_rank_1();
		test_fwd_search();
	}

	if(showSize) cout << " *** RMM total size: " << sizeRMM << " Bytes = " << (float)sizeRMM/(1024.0*1024.0) << " MB" << endl;
}

RangeMMTree64::RangeMMTree64(ulong *PSeq, ulong len, uchar *labelSeq, bool showSize) {
	// size for all tables: T_MIN_FWDI[0..255] + T_SUM_BLOCKI[0..255] + T_FWD_D[0..255][0..7] + prev_tab[0..255]
	sizeRMM =  2*256*sizeof(short int) + 256*8*sizeof(uchar) + 256*sizeof(uchar);

	// size for all variables...
	sizeRMM += 12*sizeof(ulong) + 5*sizeof(uint) + sizeof(bool);
	if (showSize) cout << " ** size for all tables and variables " << sizeRMM << endl;

	ulong nb;
	P = PSeq;
	labels = labelSeq;
	n = len;

	nb = n / S;
	if (nb%2 && nb > 1)
		nb--;
	nBin = nb*S;
	nW = n/W64;
	if (n%W64)
		nW++;
	if (TRACE) cout << " Create RangeMMTree64 with length N = 2n = " << n << ", nBin: " << nBin << ", nW: " << nW << endl;
	if (nBin >= n){
		nBin = n;
		lenLB = 0;
	}else
		lenLB = n - nBin;

	this->leaves = nBin/S;
	this->h = ceilingLog64(this->leaves, 2);
	this->firstLeaf = (ulong)pow(2, (double)h) - 1;	// the left-most leaf

	if (TRACE){
		ulong i;
		cout << "___________________________________________________________" << endl;
		cout << "P_bin :";
		for (i=0; i<nBin; i++){
			cout << readBit64(P, i);
			if ((i+1)%S == 0)
				cout << "-";
		}
		cout << "...-";
		for (; i<n; i++)
			cout << readBit64(P, i);
		cout << endl;
	}

	ulong sizeDS = n>>BW64; // size of topology
	if (n%W64)
		sizeDS++;
	sizeRMM += sizeDS*sizeof(ulong);
	if (showSize) cout << " ** size of topology " << sizeDS*sizeof(ulong) << endl;

	createMinMaxTree(showSize);

	if (TRACE){
		printTree();
		cout << "# blocks " << nb << endl;
		cout << "rank1_Bin (1's) " << rank1_Bin << endl;
		cout << "last bit of binary tree: " << nBin << endl;
		cout << "Length of the last block " << lenLB << endl;
		cout << "Sum relative to binary tree: " << sumAtPos(nBin) << endl;
		cout << "======================"<< endl;
	}

	if (RUNTEST){
		cout << "running test... " << endl;
		test_sumAtPos();
		test_select_0();
		test_selectNext0();
		test_rank_1();
		test_fwd_search();
	}

	if (showSize) cout << "*** RMM total size: " << sizeRMM << " Bytes = " << (float)sizeRMM/(1024.0*1024.0) << " MB" << endl;
}

void RangeMMTree64::createMinMaxTree(bool showSize){
	ulong groups, leavesUp, segment;
	ulong miniFwd, i, j, pos, father, node, cont, child, rb;
	long int currentSum;
	ulong *auxL, *auxR;

	groups = (uint)pow(2, (double)h-1.0);	// number of nodes at level (h-1)
	leavesBottom = 2*(leaves - groups);
	leavesUp = leaves - leavesBottom;
	firstLeaf = (uint)pow(2, (double)h) - 1;
	cantIN = firstLeaf - leavesUp;
	cantN = cantIN + leaves;

	auxL = new ulong[cantIN+leaves];	// these are auxiliary vectors to facility the compute of intervals in each internal nodes. These will be delete later.
	auxR = new ulong[cantIN+leaves];

	if (TRACE){
		cout << "leaves: " << leaves << endl;
		cout << "leaves Up: " << leavesUp << endl;
		cout << "leaves Bottom: " << leavesBottom << endl;
		cout << "internal nodes: " << cantIN << endl;
		cout << "first leaf: " << firstLeaf << endl;
		cout << "total nodes: " << cantN << endl;
	}

	// Step 1: process the n bits in groups of size s, and create the leaves in each array:	intervals of relatives excess per block, in MinLeaf and MaxLeaf,
	//         and the relative sum of excess per block, in ELeaf
	// 'i' is for bits in P, 'cont' is for count the size of each block, and 'node' is the current leaf position to set min-max values...
	j = firstLeaf; // this is for auxiliary vectors
	for (i=node=cont=0; i<nBin; i++){
		if (i== leavesBottom*S){
			node = leavesBottom;	// we move up one level
			j = cantIN;
		}
		if(cont==0){
			auxL[j] = i;
			pos = j;
			// set left boundaries, when the index j is the first child.
			while(pos && (pos+1)%2==0){
				father = (pos+1)/2;
				if ((pos+1)%2 > 1)
					father++;
				father--;
				auxL[father] = i;
				pos = father;
			}
			cont=1;
		}else{
			cont++;
			if(cont == S){
				auxR[j] = i;
				pos = j;
				// set right boundaries, when the index j is the last child. The last child always has rest == 1.
				while(pos && (pos+1)%2==1){
					father = (pos+1)/2;
					if ((pos+1)%2 > 1)
						father++;
					father--;
					auxR[father] = i;
					pos = father;
				}
				node++;
				j++;
				cont = 0;
			}
		}
	}
	if (cont)
		auxR[node] = i;

	// Step 3: create arrays of relative and super blocks...
	createTables(showSize);
	if (TRACE){
		cout << endl << "min-max Intervals..." << endl;
		for (i=0; i<(cantIN+leaves); i++)
			cout << "[" << auxL[i] << "," << auxR[i] << "] ";
		cout << endl;
	}

	ulong *Aux_Fwd_MinIN = new ulong[cantIN];
	ulong MIN_FWD;		// the lowest value for forward interval
	MIN_FWD = 0;

	// Step 3: set the min - max relatives values to each internal node.
	for(i=cantIN; i>0; i--){
		child = 2*i;		// child is the second child of node i
		if (child > cantIN){
			// here we are in a leaf
			if (child > firstLeaf)
				node = child - firstLeaf;
			else
				node = child - cantIN + leavesBottom;

			currentSum = miniFwd = 0;
			pos=(node-1)<<PotS;
			rb=(pos%W64)/BitB;
			for (j=0; j<N8S; j++){
				segment = (P[(pos+BitB*j)>>BW64] & RMM_Masks[rb]) >> (W64m8-BitB*rb);
				if (currentSum + T_MIN_FWDI[segment] > (long int)miniFwd)
					miniFwd = currentSum + T_MIN_FWDI[segment];
				currentSum += T_SUM_BLOCKI[segment];
				if (rb == N8W64-1) rb=0;
				else rb++;
			}
			pos+=S;
			for (j=0; j<N8S; j++){
				segment = (P[(pos+BitB*j)>>BW64] & RMM_Masks[rb]) >> (W64m8-BitB*rb);
				if (currentSum + T_MIN_FWDI[segment] > (long int)miniFwd)
					miniFwd = currentSum + T_MIN_FWDI[segment];
				currentSum += T_SUM_BLOCKI[segment];
				if (rb == N8W64-1) rb=0;
				else rb++;
			}
		}else{
			currentSum = sumAtPos(auxR[child-1]);
			if(auxL[child-1])
				currentSum -= sumAtPos(auxL[child-1]-1);

			miniFwd = Aux_Fwd_MinIN[child-1];
			if ((long int)Aux_Fwd_MinIN[child]-currentSum > (long int)miniFwd)
				miniFwd = Aux_Fwd_MinIN[child]-currentSum;
		}

		Aux_Fwd_MinIN[i-1] = miniFwd;
		if (miniFwd > MIN_FWD)
			MIN_FWD = miniFwd;
	}
	delete [] auxL;
	delete [] auxR;

	ulong cMIN_FWD, cMAX_BCK;
	cMIN_FWD = cMAX_BCK = 0;

	lgMIN_FWD = ceilingLog64(MIN_FWD+1, 2);
	if (TRACE)
		cout << "MIN_FWD = " << MIN_FWD << ", lgMIN_FWD = " << lgMIN_FWD << endl;

	ulong sizeMin = cantIN*lgMIN_FWD/W64;
	if ((cantIN*lgMIN_FWD)%W64)
		sizeMin++;
	Fwd_MinIN = new ulong[sizeMin];
	sizeMin *= sizeof(ulong);
	sizeRMM += sizeMin;
	if (showSize)
		cout << " ** size of Fwd_MinIN " << sizeMin << endl;

	for(i=0; i<cantIN; i++){
		setNum64(Fwd_MinIN, cMIN_FWD, lgMIN_FWD, Aux_Fwd_MinIN[i]);
		cMIN_FWD += lgMIN_FWD;
	}
}


void RangeMMTree64::createTables(bool showSize){
	ulong sizeDS;
	ulong i, j, jS, jR, cont, rb, segment;
	long int sum, sumBlock;
	long int minRelative, maxRelative;
	lenSB = leaves/2 + 1;
	long int *AuxTSBlock = new long int[lenSB];
	long int *AuxTRBlock = new long int[lenSB];

	cont = sum = minRelative = maxRelative = 0;
	jS = jR = 1;
	AuxTSBlock[0] = AuxTRBlock[0] = MAX_RelB = MAX_SupB = 0;
	for (i=rb=0; i<leaves; i++){
		sumBlock = 0;
		for (j=0; j<N8S; j++){
			segment = (P[(i*S+BitB*j)/W64] & RMM_Masks[rb]) >> (W64m8-BitB*rb);
			//printBitsUint(segment);cout<<endl;
			sumBlock -= T_SUM_BLOCKI[segment];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}
		sum += sumBlock;
		cont++;

		if (cont%2){
			// relative block
			AuxTRBlock[jR] = sumBlock;
			if (sumBlock < minRelative)
				minRelative = sumBlock;
			if (sumBlock > maxRelative)
				maxRelative = sumBlock;
			jR++;
		}else{
			// super block
			AuxTSBlock[jS] = sum;
			if (sum > (long int)MAX_SupB)
				MAX_SupB = sum;
			jS++;
			cont = 0;
		}
	}
	rank1_Bin = (nBin-sum)/2 + sum;

	if (minRelative*-1 > maxRelative)          // TAMBIEN VALIDAD EN RMM_RMQ !!
		MAX_RelB = minRelative*-1;
	else
		MAX_RelB = maxRelative;
	lgMAX_RelB = ceilingLog64(MAX_RelB+1, 2)+1;	// the first bit is for sign (1: positive; 0: negative)
	lgMAX_SupB = ceilingLog64(MAX_SupB+1, 2)+1;
	if (lgMAX_RelB <= 1)
		lgMAX_RelB = 2;
	if (lgMAX_SupB <= 1)
		lgMAX_SupB = 2;

	cont = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		cont++;
	TSBlock = new ulong[cont];
	sizeDS = cont*sizeof(ulong);
	sizeRMM += sizeDS;
	if (showSize) cout << " ** size of TSBlock[] " <<  sizeDS << " Bytes" << endl;
	for (i=cont=0; i<lenSB; i++, cont+=lgMAX_SupB){
		if(AuxTSBlock[i] >= 0)
			setBit64(TSBlock, cont);
		else{
			cleanBit64(TSBlock, cont);
			(AuxTSBlock[i]) *= -1;
		}
		setNum64(TSBlock, cont+1, lgMAX_SupB-1, AuxTSBlock[i]);
	}
	if (lenSB)
		delete [] AuxTSBlock;

	cont = lenSB*lgMAX_RelB/W64;
	if ((lenSB*lgMAX_RelB)%W64)
		cont++;
	TRBlock = new ulong[cont];
	sizeDS = cont*sizeof(ulong);
	sizeRMM += sizeDS;
	if (showSize) cout << " ** size of TRBlock[] " <<  sizeDS << " Bytes" << endl;
	for (i=cont=0; i<lenSB; i++, cont+=lgMAX_RelB){
		if(AuxTRBlock[i] >= 0)
			setBit64(TRBlock, cont);
		else{
			cleanBit64(TRBlock, cont);
			(AuxTRBlock[i]) *= -1;
		}
		setNum64(TRBlock, cont+1, lgMAX_RelB-1, AuxTRBlock[i]);
	}

	if (lenSB)
		delete [] AuxTRBlock;

	if (TRACE){
		cout << "MAX_RelB " << MAX_RelB << ", lgMAX_RelB " << lgMAX_RelB <<endl;
		cout << "MAX_SupB " << MAX_SupB << ", lgMAX_SupB " << lgMAX_SupB <<endl;
		cout << "TSBlock[1.." <<lenSB<< "]..." << endl;
		for (i=0; i<lenSB; i++){
			if (readBit64(TSBlock, i*lgMAX_SupB)==false)
				cout << "-";
			cout << getNum64(TSBlock, i*lgMAX_SupB+1, lgMAX_SupB-1) << " ";
		}
		cout << endl;
		cout << "TRBlock[1.." <<lenSB<< "]..." << endl;
		for (i=0; i<lenSB; i++){
			if (readBit64(TRBlock, i*lgMAX_RelB)==false)
				cout << "-";
			cout << getNum64(TRBlock, i*lgMAX_RelB+1, lgMAX_RelB-1) << " ";
		}
		cout << endl;
	}
}

// here, 1 gives a positive increment and 0 gives a negative increment
long int RangeMMTree64::sumAtPos(ulong pos){
	ulong rb, q, i, l, k, j=pos+1;
	ulong blk = j>>PotS;
	ulong srBlock = blk>>1;
	long int aux, sum = getNum64(TSBlock, srBlock*lgMAX_SupB+1, lgMAX_SupB-1);
	if (!readBit64(TSBlock, srBlock*lgMAX_SupB))
		sum *= -1;

	if(blk%2){
		srBlock++;
		if (srBlock<lenSB){
			aux = getNum64(TRBlock, srBlock*lgMAX_RelB+1, lgMAX_RelB-1);
			if (!readBit64(TRBlock, srBlock*lgMAX_RelB))
				aux *= -1;
			sum += aux;
		}else
			blk--;
	}

	if (j>=BitB){
		l=blk<<PotS;
		rb=(l%W64)/BitB;
		for (k=j-BitB; l<=k; l+=BitB){
			q = (P[l>>BW64] & RMM_Masks[rb]) >> (W64m8-BitB*rb);
			sum -= T_SUM_BLOCKI[q];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}

	}
	if (j%BitB){
		for (i=pos-j%BitB+1; i<j; i++){
			if (readBit64(P, i))
				sum++;
			else
				sum--;
		}
	}
	return sum;
}

void RangeMMTree64::test_sumAtPos(){
	ulong i, j, k;
	long int sum, sumPos;

	cout << "RangeMMTree64::test_sumAtPos..." << endl;
	/*i=160;
	for (j=sum=0; j<=i; j++){
		if(readBit64(P, j))
			sum++;
		else
			sum--;
	}
	sumPos = sumAtPos(i);
	if (sum != sumPos){
		cout << "ERROR !! sumAtPos(" << i << ") = " << sumPos << " != sum = " << sum << endl;
		exit(1);
	}
	exit(0);*/

	for (k=0; k<TEST; k++){
		i = (rand() % (n-2));
		for (j=sum=0; j<=i; j++){
			if(readBit64(P, j))
				sum++;
			else
				sum--;
		}
		sumPos = sumAtPos(i);
		if (sum != sumPos){
			cout << "ERROR !! sumAtPos(" << i << ") = " << sumPos << " != sum = " << sum << endl;
			exit(1);
		}
	}
	cout << "  test_sumAtPos OK !!" << endl;
}

void RangeMMTree64::printTree(){
	ulong i, j, cant, acum, rb, segment;
	long int currSum, mini;
	ulong cMIN_FWD, cMIN_BCK, cMAX_BCK, cMAX_FWD;
	cMIN_FWD = cMIN_BCK = cMAX_BCK = cMAX_FWD = 0;

	cout << "Min-max Binary Tree...i(sum)[minFwd]" << endl;
	cant = 1;
	cout << "Fwd_MinIN[] = ";
	for (acum=j=i=0; i<cantIN; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " - ";
			acum = 0;
		}

		mini = getNum64(Fwd_MinIN, cMIN_FWD, lgMIN_FWD);
		cout << i << "_IN_Fwd[" << mini << "] ";
		cMIN_FWD += lgMIN_FWD;
	}
	cout << endl;

	i=leavesBottom;
	rb=((i<<PotS)%W64)/BitB;
	for (; i<leaves; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " - ";
			acum = 0;
		}
		currSum = 0;
		mini = 0;
		for (j=0; j<N8S; j++){
			segment = (P[(i*S+BitB*j)>>BW64] & RMM_Masks[rb]) >> (W64m8-BitB*rb);
			if (currSum + T_MIN_FWDI[segment] > mini)
				mini = currSum + T_MIN_FWDI[segment];
			currSum += T_SUM_BLOCKI[segment];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}
		cout << i << "_h_Fwd[" << mini << "](" << currSum << ")";
	}
	cout << endl;

	i=0;
	rb=((i<<PotS)%W64)/BitB;
	for (; i<leavesBottom; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " - ";
			acum = 0;
		}

		currSum = 0;
		mini = 0;
		for (j=0; j<N8S; j++){
			segment = (P[(i*S+BitB*j)>>BW64] & RMM_Masks[rb]) >> (W64m8-BitB*rb);
			if (currSum + T_MIN_FWDI[segment] > mini)
				mini = currSum + T_MIN_FWDI[segment];
			currSum += T_SUM_BLOCKI[segment];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}
		cout << i << "_h_Fwd[" << mini << "](" << currSum << ")";
	}
	cout << endl;
}

// for 0 <= i < n
ulong RangeMMTree64::binRank_1(ulong i){
	long int sum = sumAtPos(i);
	return sum+((i+1-sum)>>1);
}

// for 0 <= i < n
ulong RangeMMTree64::binRank_0(ulong i){
	long int sum = sumAtPos(i);
	return ((i-sum)>>1)-sum;
}

// for 0 <= i < n
ulong RangeMMTree64::rank_1(ulong i){
	if(i < nBin)
		return binRank_1(i);
	else{
		if (nBin == i)
			return rank1_Bin + readBit64(P,i);

		if(i > n-1)
			i=n-1;

		ulong rank=rank1_Bin;
		ulong b, rest, q, x = nBin;

		while(x+BmOne <= i){
			b = x>>BW64;
			rest = (x+BitB)%W64;
			q = (P[b] >> (W64-rest)) & 0xff;
			//printBitsUint(q);cout<<endl;
			rank += __popcount_tab[q];
			x += BitB;
		}

		// check last segment (< S) bit by bit...
		while (x<=i){
			if(readBit64(P,x)) rank++;
			x++;
		}

		return rank;
	}
}

void RangeMMTree64::test_rank_1(){
	ulong sum,rank,i,j,k;

	cout << "RangeMinMaxTree_Bin::test_rank_1..." << endl;
	/*i=58;
	for (j=sum=0; j<=i; j++){
		if(readBit64(P, j))
			sum++;
	}
	cout << "rank_1("<<i<<") = " << sum << endl;
	rank = rank_1(i);
	cout << "function rank_1("<<i<<") = " << rank << endl;
	exit(0);*/

	for (k=0; k<TEST; k++){
		i = (rand() % (n-10));
		for (j=sum=0; j<=i; j++){
			if(readBit64(P, j))
				sum++;
		}
		rank = rank_1(i);
		if (sum != rank){
			cout << "ERROR !! rank1(" << i << ") = " << rank << " != sum = " << sum << endl;
			exit(1);
		}
	}
	cout << "  test_rank_1 OK !!" << endl;
}


ulong RangeMMTree64::binSelect_0(ulong i){
	ulong l, r, m, bit, prev, sbm, q, rest, ze;
	long int aux;

	m = lenSB>>1;
	if(m){
		l=1; r=lenSB-1;
		bit = (m<<1)<<PotS;
		aux = getNum64(TSBlock, m*lgMAX_SupB+1, lgMAX_SupB-1);
		if (!readBit64(TSBlock, m*lgMAX_SupB))
			aux *= -1;
		sbm = (bit-aux)>>1;
		if (m>1){
			aux = getNum64(TSBlock, (m-1)*lgMAX_SupB+1, lgMAX_SupB-1);
			if (!readBit64(TSBlock, (m-1)*lgMAX_SupB))
				aux *= -1;
			prev = (bit-(S<<1)-aux)>>1;
		}else
			prev=0;

		while(i<=prev || i>sbm){
			if (i > sbm)
				l=m+1;
			else
				r=m-1;

			if (l<=r){
				m = l+((r-l)>>1);
				bit = (m<<1)<<PotS;
				aux = getNum64(TSBlock, m*lgMAX_SupB+1, lgMAX_SupB-1);
				if (!readBit64(TSBlock, m*lgMAX_SupB))
					aux *= -1;
				sbm = (bit-aux)>>1;
				if (m>1){
					aux = getNum64(TSBlock, (m-1)*lgMAX_SupB+1, lgMAX_SupB-1);
					if (!readBit64(TSBlock, (m-1)*lgMAX_SupB))
						aux *= -1;
					prev = (bit-(S<<1)-aux)>>1;
				}else
					prev = 0;
			}else break;
		}

		// m is the super block where i-th zero is
		aux = getNum64(TRBlock, m*lgMAX_RelB+1, lgMAX_RelB-1);
		if (!readBit64(TRBlock, m*lgMAX_RelB))
			aux *= -1;
		l = (S-aux)>>1;
		if (i > prev+l){
			bit -= S;
			sbm = prev+l;
		}else{
			bit -= S<<1;
			sbm = prev;
		}
	}else{
		bit = 0;
		sbm = 0;
	}

	ulong b = bit>>BW64;
	rest = (bit+BitB)%W64;
	q = (P[b] >> (W64-rest)) & 0xFFUL;
	//printBitsUint(P[b]);cout<<endl;
	//printBitsUint(q);cout<<endl;
	ze = BitB-popcount_Rank64(q);

	while(sbm+ze < i){
		sbm += ze;
		bit += BitB;
		b = bit>>BW64;
		rest = (bit+BitB)%W64;
		q = (P[b] >> (W64-rest)) & 0xFFUL;
		//printBitsUint(q);cout<<endl;
		ze = BitB-popcount_Rank64(q);
	}

	// check last segment (< S) bit by bit...
	while (sbm < i){
		if(!readBit64(P,bit)) sbm++;
		bit++;
	}

	return bit-1;
}

// for 0 <= i < n/2
ulong RangeMMTree64::select_0(ulong i){
	if(i <= (n>>1)){
		ulong ze = nBin - rank1_Bin;

		if(i <= ze)
			return binSelect_0(i);
		else{
			ulong sel=ze;
			ulong b, rest, q, x = nBin;

			b = x>>BW64;
			rest = BitB;
			q = (P[b] >> (W64-rest)) & 0xFFUL;
			ze = BitB-popcount_Rank64(q);

			while(sel+ze < i){
				sel += ze;
				x += 8;
				b = x>>BW64;
				rest = (x+BitB)%W64;
				q = (P[b] >> (W64-rest)) & 0xFFUL;
				//printBitsUint(q);cout<<endl;
				ze = BitB-popcount_Rank64(q);
			}

			// check last segment (< S) bit by bit...
			while (sel < i){
				if(!readBit64(P,x)) sel++;
				x++;
			}

			return x-1;
		}
	}
	return 0;
}

void RangeMMTree64::test_select_0(){
	ulong sum0,sel,i,j,k;

	cout << "RangeMMTree64::test_select_0..." << endl;
	/*i=29;
	for (j=sum0=0; sum0<i; j++){
		if(!readBit64(P, j))
			sum0++;
	}
	j--;
	cout << "brute select_0("<<i<<") = " << j << endl;
	sel = select_0(i);
	cout << "funt. select_0 = " << sel << endl;
	exit(0);*/

	for (k=0; k<TEST; k++){
		i = (rand() % (n/2));
		for (j=sum0=0; sum0<i; j++){
			if(!readBit64(P, j))
				sum0++;
		}
		j--;
		sel = select_0(i);
		if (j != sel){
			cout << "ERROR !! sel0(" << i << ") = " << sel << " != select_brute = " << j << endl;
			exit(1);
		}
	}
	cout << "  test_select_0 OK !!" << endl;
}



void RangeMMTree64::test_selectNext0(){
	long int sel;
	ulong i,j,k;

	cout << "RangeMMTree64::test_next_0..." << endl;
	//sel = selectNext0(92);
	//cout << "next_0(92) = " << sel << endl;

	for (k=0; k<TEST; k++){
		j = (rand() % (n-1));
		for (i=j; i<n && readBit64(P, i); i++);
		sel = selectNext0(j);
		if (sel != (long int)i){
			cout << "ERROR !! next_0(" << j << ") = " << sel << " != " << i << endl;
			exit(1);
		}

	}
	cout << "  test_next_0 OK !!" << endl;
}

ulong RangeMMTree64::selectNext0(const ulong count) const
{
	ulong des,aux2;
	des=count%W64;
	aux2= ~(P[count>>BW64]) << des;
	if (aux2 > 0) {
		     if ((aux2&0xff00000000000000UL) > 0) return count+8-prev_tab[aux2>>56];
		else if ((aux2&0xff000000000000UL) > 0) return count+16-prev_tab[aux2>>48];
		else if ((aux2&0xff0000000000UL) > 0) return count+24-prev_tab[aux2>>40];
		else if ((aux2&0xff00000000UL) > 0) return count+32-prev_tab[aux2>>32];
		else if ((aux2&0xff000000UL) > 0) return count+40-prev_tab[aux2>>24];
		else if ((aux2&0xff0000UL) > 0) return count+48-prev_tab[aux2>>16];
		else if ((aux2&0xff00UL) > 0) return count+56-prev_tab[aux2>>8];
		else {return count+64-prev_tab[aux2];}
	}

	for (ulong i=(count>>BW64)+1;i<nW;i++) {
		aux2=~(P[i]);
		if (aux2 > 0) {
			     if ((aux2&0xff00000000000000UL) > 0) return (i<<BW64)+8-prev_tab[aux2>>56];
			else if ((aux2&0xff000000000000UL) > 0) return (i<<BW64)+16-prev_tab[aux2>>48];
			else if ((aux2&0xff0000000000UL) > 0) return (i<<BW64)+24-prev_tab[aux2>>40];
			else if ((aux2&0xff00000000UL) > 0) return (i<<BW64)+32-prev_tab[aux2>>32];
			else if ((aux2&0xff000000UL) > 0) return (i<<BW64)+40-prev_tab[aux2>>24];
			else if ((aux2&0xff0000UL) > 0) return (i<<BW64)+48-prev_tab[aux2>>16];
			else if ((aux2&0xff00UL) > 0) return (i<<BW64)+56-prev_tab[aux2>>8];
			else {return (i<<BW64)+64-prev_tab[aux2];}
		}
	}
	return n;
}


// return true and the position 'pos', pos < ini+len, where the sum from ini to ini+len is d. If not found then return false
bool RangeMMTree64::fwd_block(ulong x, ulong *d, ulong* pos){
	ulong rest, b, q, len = S - x%S;
	ulong tar = *d;

	while(len > BmOne){
		b = x>>BW64;
		rest = (x+BitB)%W64;
		if (b == (x+BmOne)>>BW64)			// x and (x+S) are in the same word...
			q = (P[b] >> (W64-rest)) & 0xFFUL;
		else
			q = ((P[b] << rest) | (P[b+1] >> (W64-rest))) & 0xFFUL;
		//printBitsUint(q);cout<<endl;

		if (tar > T_MIN_FWDI[q])
			tar -= T_SUM_BLOCKI[q];
		else{
			// here is the target...
			*pos = x+T_FWD_D[q][T_MIN_FWDI[q]-tar];
			return true;
		}
		x += BitB;
		len -= BitB;
	}
	if (len){
		// check last segment (len < 8) bit by bit...
		long int sum;
		if (readBit64(P,x)) sum = -1;
		else sum = 1;
		while (sum!=(long int)tar && len>1){
			x++;
			if(readBit64(P,x)) sum--;
			else sum++;
			len--;
		}
		if (sum == (long int)tar){
			*pos = x;
			return true;
		}
		*d = tar-sum;
	}else
		*d = tar;

	return false;
}

// search in the last segment of length (n-x)
bool RangeMMTree64::fwd_Lastblock(ulong x, ulong d, ulong *pos){
	ulong b, rest, q;

	while(x+BmOne < n){
		b = x>>BW64;
		rest = (x+BitB)%W64;
		if (b == (x+BmOne)>>BW64)			// x and (x+S) are in the same word...
			q = (P[b] >> (W64-rest)) & 0xFFUL;
		else
			q = ((P[b] << rest) | (P[b+1] >> (W64-rest))) & 0xFFUL;
		//printBitsUint(q);cout<<endl;

		if (d > T_MIN_FWDI[q])
			d -= T_SUM_BLOCKI[q];
		else{
			// here is the target...
			*pos = x+T_FWD_D[q][T_MIN_FWDI[q]-d];
			return true;
		}
		x += BitB;
	}

	// check last segment (< S) bit by bit...
	//printBitsUint(P[x/W]);cout<<endl;
	long int sum;
	if (readBit64(P,x)) sum = -1;
	else sum = 1;
	while (sum!=(long int)d && x<n){
		x++;
		if(readBit64(P,x)) sum--;
		else sum++;
	}
	if (x < n && sum == (long int)d){
		*pos = x;
		return true;
	}

	return false;
}

// return true and the position pos >= ini that the sum from ini to pos is d. If not found then return false
//d  Always is grater than 0
bool RangeMMTree64::binFwd_search(ulong ini, ulong *d, ulong* pos){
	ulong node = ini>>PotS;				// block for 'ini'
	ulong mini, j, rb, q, height, dist = 1;

	// [1]- Check if the answer is inside the block...
	if (fwd_block(ini, d, pos))
		return true;
	else{
		if (leaves == node+1)
			return false;
	}

	// [2]- We looking for in the brother leaf (at the same super block) this target 'target'
	ulong target = *d;
	if (node%2 == 0){
		node++;
		ini = node<<PotS;
		rb=(ini%W64)/BitB;
		for (j=0; j<N8S; j++){
			q = (P[(ini+BitB*j)>>BW64] & RMM_Masks[rb]) >> (W64m8-BitB*rb);
			if (target > T_MIN_FWDI[q])
				target -= T_SUM_BLOCKI[q];
			else{
				// here is the target...
				*pos = ini+BitB*j+T_FWD_D[q][T_MIN_FWDI[q]-target];
				return true;
			}
			if (rb == N8W64-1) rb=0;
			else rb++;
		}
	}
	if (node == leaves-1){
		*d = target;
		return false;
	}

	// [3]- We climb recomputing the target and looking for the node that contains this target.
	height = h-2;
	if (node < leavesBottom){
		if (node == leavesBottom-1)
			node = cantIN-1;
		else{
			node += firstLeaf;
			height++;
		}
	}else
		node += cantIN-leavesBottom;
	node >>= 1;
	mini = getNum64(Fwd_MinIN, node*lgMIN_FWD, lgMIN_FWD);
	while (target > mini){
		target -= sumOfNode(node);
		if (node%2)
			node++;
		else{
			float lg2 = log(node+2)/log(2);
			if(node+2 == pow(2,lg2)){ // this line check id the node is a rightmost node for its level
				*d = target;
				return false;
			}
			node>>=1;				// go to my uncle
			height--;
			dist++;
		}
		mini = getNum64(Fwd_MinIN, node*lgMIN_FWD, lgMIN_FWD);
	}

	// [4]- We move down recomputing the target and reach to the correct leaf that has the target found.
	node = (node<<1) + 1;
	dist--;
	while (node < cantIN){
		mini = getNum64(Fwd_MinIN, node*lgMIN_FWD, lgMIN_FWD);
		if (target > mini){
			target -= sumOfNode(node);
			node++;
		}
		node = (node<<1) + 1;
		dist--;
	}
	if (node > firstLeaf)						// at this point 'node' is a leaf, node >= firstLeaf
		node -= firstLeaf;
	else
		node += leavesBottom - cantIN;

	ini = node<<PotS;
	rb=(ini%W64)/BitB;
	for (j=0; j<N8S; j++){
		q = (P[(ini+BitB*j)>>BW64] & RMM_Masks[rb]) >> (W64m8-BitB*rb);
		if (target > T_MIN_FWDI[q])
			target -= T_SUM_BLOCKI[q];
		else{
			// here is the target...
			*pos = ini+BitB*j+T_FWD_D[q][T_MIN_FWDI[q]-target];
			return true;
		}
		if (rb == N8W64-1) rb=0;
		else rb++;
	}

	// [5]- We return the answer inside the leaf that contains the target
	return fwd_block(ini+S, &target, pos);
}

bool RangeMMTree64::fwd_search(ulong x, ulong *d, ulong *pos){
	if(x < nBin){
		if (binFwd_search(x, d, pos))
			return true;
		else
			return fwd_Lastblock(nBin, *d, pos);

	}else
		return fwd_Lastblock(x, *d, pos);
}

void RangeMMTree64::test_fwd_search(){
	ulong i, j, r, pos, d, tar;
	long int sum;

	cout << "RangeMMTree64::test_fwd_search..." << endl;

	/*d=1;
	r=5;
	for (sum=0; r<n; r++){
		if(readBit64(P, r))
			sum--;
		else
			sum++;
		if (sum == (long int)d)
			break;
	}
	j=5;
	cout << "fwd_search(" << j << ", " << d << ") = " << r << endl;
	fwd_search(j, &d, &pos);
	if (pos != r)
		cout << "ERROR !! pos = " << pos << " != r = " << r << endl;
	exit(0);*/

	for (d=1; d<=15; d++){
		for (i=0; i<TEST; i++){
			j = (rand() % (n-1));

			for (sum=0, r=j; r<n; r++){
				if(readBit64(P, r))
					sum--;
				else
					sum++;
				if (sum == (long int)d)
					break;
			}
			tar = d;
			if (fwd_search(j, &tar, &pos)){
				if (r==n){
					cout << "fwd_search(" << j << ", " << d << ")" << endl;
					cout << "ERROR !! pos = " << pos << " and NOT FOUND" << endl;
					exit(1);
				}
				if (pos != r){
					cout << "fwd_search(" << j << ", " << d << ")" << endl;
					cout << "ERROR !! pos = " << pos << " != r = " << r << endl;
					exit(1);
				}
			}else{
				if (r < n){
					cout << "fwd_search(" << j << ", " << d << ")" << endl;
					cout << "ERROR !! fwd_search is FALSE and FOUND in r = " << r << endl;
					exit(1);
				}
			}
		}
	}
	cout << "  test_fwd_search OK !!" << endl;
}

// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
long int RangeMMTree64::sumOfNode(ulong node){
	uint ini, end;

	ini = end = node;
	while(ini<cantIN)
		ini=(ini<<1)+1;
	if (ini >= firstLeaf)
		ini -= firstLeaf;
	else
		ini = ini-cantIN+leavesBottom;

	while(end<cantIN)
		end=(end+1)<<1;
	if(end >= firstLeaf)
		end -= firstLeaf-1;
	else
		end = end-cantIN+leavesBottom+1;

	end>>=1;
	long int sumE = getNum64(TSBlock, end*lgMAX_SupB+1, lgMAX_SupB-1);
	if (!readBit64(TSBlock, end*lgMAX_SupB))
		sumE *= -1;
	if (ini==0)
		return sumE;

	ini>>=1;
	long int sumI = getNum64(TSBlock, ini*lgMAX_SupB+1, lgMAX_SupB-1);
	if (!readBit64(TSBlock, ini*lgMAX_SupB))
		sumI *= -1;

	return (sumI-sumE);
}


// return true and the relative position of the child (in 'child') labeled with 'c' for the node 'x'. The parameter 'end' is the position of the last child of x
bool RangeMMTree64::thereIsChild(ulong x, uchar c, ulong *child, uint len){
	uint m;
	ulong ini = rank_1(x-1) - 1;
	ulong end = ini+len-1;
	*child = ini;
	// binary search in the interval LbRev[ini..end]

	while(ini <= end){
		m = ini+((end-ini)>>1);
		if (labels[m] == c){
			*child = m - *child + 1;
			return true;
		}
		if (labels[m] < c)
			ini=m+1;
		else{
			if(m)
				end=m-1;
			else
				return false;
		}
	}
	if (ini <= end && labels[ini] == c){
		*child = ini - *child + 1;
		return true;
	}
	return false;
}

// return true if this fictitious node has a unique child
bool RangeMMTree64::hasUniqueChild(ulong x){
	if (selectNext0(x) == x+1)
		return true;
	return false;
}

// it gives the size of the subtree rooted at x (without fictitious nodes)
ulong RangeMMTree64::subTreeSize(ulong x){
	if (readBit64(P,x)){
		ulong pos, len=0;
		ulong d=1;
		if(fwd_search(x, &d, &pos))
			len = ((pos-x)>>1) + 1;

		return len;
	}else
		return 1;	// this node is a leaf
}

// save the Data Structure in file 'fileName' and return the amount of bytes saved
ulong RangeMMTree64::saveDS(char *fileName, bool showSize){
	cout << "Save data structure in " << fileName << endl;
	ofstream os (fileName, ios::binary);
	cout << "   Saving. Data structure size: " << sizeRMM << endl;

	if(TRACE){
		cout << "Variables load... " << endl;
		cout << "n " << n << endl;
		cout << "nW " << nW << endl;
		cout << "rank1_Bin " << rank1_Bin << endl;
		cout << "nBin " << nBin << endl;
		cout << "lenLB " << lenLB << endl;
		cout << "h " << h << endl;
		cout << "cantN " << cantN << endl;
		cout << "cantIN " << cantIN << endl;
		cout << "leaves " << leaves << endl;
		cout << "leavesBottom " << leavesBottom << endl;
		cout << "firstLeaf " << firstLeaf << endl;
		cout << "lenSB " << lenSB << endl;
		cout << "lgMIN_FWD " << lgMIN_FWD << endl;
		cout << "MAX_RelB " << MAX_RelB << endl;
		cout << "lgMAX_RelB " << lgMAX_RelB << endl;
		cout << "MAX_SupB " << MAX_SupB << endl;
		cout << "lgMAX_SupB " << lgMAX_SupB << endl;
	}

	os.write((const char*)&n, sizeof(ulong));
	os.write((const char*)&nW, sizeof(ulong));
	os.write((const char*)&rank1_Bin, sizeof(ulong));
	os.write((const char*)&nBin, sizeof(ulong));
	os.write((const char*)&lenLB, sizeof(uint));
	os.write((const char*)&h, sizeof(uint));
	os.write((const char*)&cantN, sizeof(ulong));
	os.write((const char*)&cantIN, sizeof(ulong));
	os.write((const char*)&leaves, sizeof(ulong));
	os.write((const char*)&leavesBottom, sizeof(ulong));
	os.write((const char*)&firstLeaf, sizeof(ulong));
	os.write((const char*)&lenSB, sizeof(ulong));
	os.write((const char*)&lgMIN_FWD, sizeof(uint));
	os.write((const char*)&MAX_RelB, sizeof(ulong));
	os.write((const char*)&lgMAX_RelB, sizeof(uint));
	os.write((const char*)&MAX_SupB, sizeof(ulong));
	os.write((const char*)&lgMAX_SupB, sizeof(uint));

	// size of tables: T_MIN_FWDI[0..255] + T_SUM_BLOCKI[0..255] + T_FWD_D[0..255][0..7] + prev_tab[0..255]
	ulong sizeDT =  2*256*sizeof(short int) + 256*8*sizeof(uchar) + 256*sizeof(uchar);
	// size for variables
	sizeDT += 12*sizeof(ulong) + 5*sizeof(uint) + sizeof(bool);
	if (showSize) cout << " ** size of all tables and variables " << sizeDT << endl;

	ulong size = n>>BW64;
	if (n % W64)
		size++;
	os.write((const char*)P, size*sizeof(ulong));				// save P[]
	sizeDT += size*sizeof(ulong);
	if(showSize) cout << " ** size of topology P[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = cantIN*lgMIN_FWD/W64;
	if ((cantIN*lgMIN_FWD)%W64)
		size++;
	os.write((const char*)Fwd_MinIN, size*sizeof(ulong));		// save Bkwd_MinIN[]
	sizeDT += size*sizeof(ulong);
	if(showSize) cout << " ** size of Fwd_MinIN[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		size++;
	os.write((const char*)TSBlock, size*sizeof(ulong));			// save TSBlock[]
	sizeDT += size*sizeof(ulong);
	if(showSize) cout << " ** size of TSBlock[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = lenSB*lgMAX_RelB/W64;
	if ((lenSB*lgMAX_RelB)%W64)
		size++;
	os.write((const char*)TRBlock, size*sizeof(ulong));			// save TRBlock[]
	sizeDT += size*sizeof(ulong);
	if(showSize) cout << " ** size of TRBlock[] " << size*sizeof(ulong) << " Bytes" << endl;

	os.close();
	if(showSize) cout << "   Total bytes saved: " << sizeDT << endl;
	return sizeDT;
}

// load the Data Structure from the file 'fileName'
void RangeMMTree64::loadDS(char *fileName, bool showSize){
	// size of tables: T_MIN_FWDI[0..255] + T_SUM_BLOCKI[0..255] + T_FWD_D[0..255][0..7] + prev_tab[0..255]
	sizeRMM =  2*256*sizeof(short int) + 256*8*sizeof(uchar) + 256*sizeof(uchar);

	// size for variables
	sizeRMM += 12*sizeof(ulong) + 5*sizeof(uint) + sizeof(bool);
	if (showSize) cout << " ** size of all tables and variables " << sizeRMM << endl;

	cout << " Load trie from " << fileName << endl;
	ifstream is(fileName, ios::binary);

	is.read((char*)&n, sizeof(ulong));
	is.read((char*)&nW, sizeof(ulong));
	is.read((char*)&rank1_Bin, sizeof(ulong));
	is.read((char*)&nBin, sizeof(ulong));
	is.read((char*)&lenLB, sizeof(uint));
	is.read((char*)&h, sizeof(uint));
	is.read((char*)&cantN, sizeof(ulong));
	is.read((char*)&cantIN, sizeof(ulong));
	is.read((char*)&leaves, sizeof(ulong));
	is.read((char*)&leavesBottom, sizeof(ulong));
	is.read((char*)&firstLeaf, sizeof(ulong));
	is.read((char*)&lenSB, sizeof(ulong));
	is.read((char*)&lgMIN_FWD, sizeof(uint));
	is.read((char*)&MAX_RelB, sizeof(ulong));
	is.read((char*)&lgMAX_RelB, sizeof(uint));
	is.read((char*)&MAX_SupB, sizeof(ulong));
	is.read((char*)&lgMAX_SupB, sizeof(uint));

	if(TRACE){
		cout << "Variables load... " << endl;
		cout << "n " << n << endl;
		cout << "nW " << nW << endl;
		cout << "rank1_Bin " << rank1_Bin << endl;
		cout << "nBin " << nBin << endl;
		cout << "lenLB " << lenLB << endl;
		cout << "h " << h << endl;
		cout << "cantN " << cantN << endl;
		cout << "cantIN " << cantIN << endl;
		cout << "leaves " << leaves << endl;
		cout << "leavesBottom " << leavesBottom << endl;
		cout << "firstLeaf " << firstLeaf << endl;
		cout << "lenSB " << lenSB << endl;
		cout << "lgMIN_FWD " << lgMIN_FWD << endl;
		cout << "MAX_RelB " << MAX_RelB << endl;
		cout << "lgMAX_RelB " << lgMAX_RelB << endl;
		cout << "MAX_SupB " << MAX_SupB << endl;
		cout << "lgMAX_SupB " << lgMAX_SupB << endl;
	}

	ulong sizeDS = n>>BW64;
	if (n % W64)
		sizeDS++;
	P = new ulong[sizeDS];
	is.read((char*)P, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(showSize) cout << " ** size of topology P[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	sizeDS = cantIN*lgMIN_FWD/W64;
	if ((cantIN*lgMIN_FWD)%W64)
		sizeDS++;
	Fwd_MinIN = new ulong[sizeDS];
	is.read((char*)Fwd_MinIN, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(showSize) cout << " ** size of Fwd_MinIN[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	sizeDS = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		sizeDS++;
	TSBlock = new ulong[sizeDS];
	is.read((char*)TSBlock, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(showSize) cout << " ** size of TSBlock[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	sizeDS = lenSB*lgMAX_RelB/W64;
	if ((lenSB*lgMAX_RelB)%W64)
		sizeDS++;
	TRBlock = new ulong[sizeDS];
	is.read((char*)TRBlock, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(showSize) cout << " ** size of TRBlock[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	is.close();
	if(TRACE) cout << " Data Structure loaded, total bytes:" << sizeRMM << endl;
}

RangeMMTree64::~RangeMMTree64() {
	n = nW = rank1_Bin = nBin = lenLB = h = cantN = cantIN =
	leaves = leavesBottom = firstLeaf = lenSB = lgMIN_FWD =
	MAX_RelB = lgMAX_RelB = MAX_SupB = lgMAX_SupB = sizeRMM = 0;

	delete [] Fwd_MinIN;
	delete [] TSBlock;
	delete [] TRBlock;
	cout << " ~ RangeMMTree64 destroyed !!" << endl;
}

