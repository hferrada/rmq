/*
 * RMQRMM64.cpp
 *
 *  Created on: 16-06-2014
 *      Author: hector
 */

#include "includes/RMQRMM64.h"

bool RMQRMM64::TRACE = false;
bool RMQRMM64::RUNTEST = false;
uint RMQRMM64::TEST = 10000;

RMQRMM64::RMQRMM64(char *fileName){
	loadDS(fileName);

	if (RUNTEST){
		test_sumAtPos();
		test_rank_1();
		test_select_1();
		if (leaves>0) test_rmqi();
	}
}

RMQRMM64::RMQRMM64(long int *A, ulong len) {
	if (TRACE){
		cout << "*** RMQ structure for A[0.." << len-1 << "] and leaf length S = " << Srmq << endl;
		cout << "Create BP sequence for 2D Min Heap..." << endl;
	}

	sizeRMM = 2*512 + 2048;											// size for T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[]
	sizeRMM += 10*sizeof(ulong) + 10*sizeof(uint) + sizeof(int);	// size for variables

	nP = (len+1)<<1;
	ulong lenP = nP >> BW64;
	if (nP % W64)
		lenP++;
	P = new ulong[lenP];
	sizeRMM += lenP*sizeof(ulong);
	if (TRACE) cout << " ** size of topology " << lenP*sizeof(ulong) << " Bytes" << endl;

	ulong pos = 2;
	ulong i;
	long int j;
	setBit64(P, 0);
	setBit64(P, 1);

	long int *fath = new long int[len];
	fath[0] = -1;
	for(i=1; i<len; i++){
		j=i-1;
		while(j>0 && A[j]>=A[i]){
			cleanBit64(P, pos);				// to close the nodes opened before me
			pos++;
			j = fath[j];
		}
		if (j==0){
			if(A[0] >= A[i]){
				cleanBit64(P, pos);			// to close the nodes opened before me
				pos++;
				fath[i] = -1;
			}else
				fath[i] = 0;
		}else
			fath[i] = j;
		setBit64(P, pos);	// to open this node
		pos++;
	}
	cleanBit64(P, pos);
	cleanBit64(P, pos+1);
	delete [] fath;

	if(RUNTEST){ // test of balanced sequence
		long int sum = 0;
		ulong r=0;

		for (; r<nP; r++){
			if(readBit64(P, r))
				sum++;
			else sum--;
		}
		if(sum != 0){
			cout << " ERROR. P[] is not a balanced sequence of parentheses !! " << endl;
			exit(0);
		}else
			cout << " P[] is a well balanced sequence of parentheses !! " << endl;
	}

	if (TRACE){
		cout << " BP sequence, P[0.." << nP-1 << "]" << endl;
		for (pos=0; pos<lenP; pos++)
			printBitsUlong(P[pos]);
		cout << endl;
	}

	ulong nb = nP / Srmq;
	if (nb%2 && nb > 1)
		nb--;
	nBin = nb*Srmq;
	nW = nP/W64;
	if (nP%W64)
		nW++;

	if (nBin >= nP){
		nBin = nP;
		lenLB = 0;
	}else
		lenLB = nP - nBin;

	if (TRACE)
		cout << "Create RangeMinMaxTree_Bin with length N = 2n = " << nP << ", nBin: " << nBin << ", nW: " << nW << endl;

	this->leaves = nBin/Srmq;
	this->h = ceilingLog64(this->leaves, 2);
	this->firstLeaf = (ulong)pow(2, (double)h) - 1;	// the left-most leaf

	if (TRACE){
		ulong i;
		cout << "___________________________________________________________" << endl;
		cout << "P_bin :";
		for (i=0; i<nBin; i++){
			cout << readBit64(P, i);
			if ((i+1)%Srmq == 0)
				cout << "-";
		}
		cout << endl << "Last Block :" << endl;
		for (; i<nP; i++)
			cout << readBit64(P, i);
		cout << endl;
	}

	createMinMaxTree();
	if (TRACE){
		if (leaves>0) printTree();
		cout << "# blocks " << nb << endl;
		cout << "# blocks of the binary tree " << nb << endl;
		cout << "rank1_Bin (1's) " << rank1_Bin << endl;
		cout << "last bit of binary tree: " << nBin << endl;
		cout << "Length of the last block " << lenLB << endl;
		cout << "Sum relative to binary tree: " << sumAtPos(nBin) << endl;
		cout << "======================"<< endl;
	}
	if (RUNTEST){
		test_sumAtPos();
		test_rank_1();
		test_select_1();
		if (leaves>0) test_rmqi();
	}
}

void RMQRMM64::createMinMaxTree(){
	ulong groups, leavesUp, sizeDS;
	ulong i, j, position, father, node, cont, child, rb;
	int miniBck;
	ulong *auxL, *auxR;
	int *Aux_Bkwd_MinIN;

	if (h>0){
		groups = (ulong)pow(2, (double)h-1.0);	// number of nodes at level (h-1)
		leavesBottom = (leaves - groups)<<1;
		leavesUp = leaves - leavesBottom;
	}else{
		leavesBottom = leaves;
		leavesUp = 0;
	}

	firstLeaf = (ulong)pow(2, (double)h) - 1;
	cantIN = firstLeaf - leavesUp;
	cantN = cantIN + leaves;

	if (TRACE && leaves>0){
		cout << "leaves: " << leaves << endl;
		cout << "leaves Up: " << leavesUp << endl;
		cout << "leaves Bottom: " << leavesBottom << endl;
		cout << "internal nodes: " << cantIN << endl;
		cout << "first leaf: " << firstLeaf << endl;
		cout << "total nodes: " << cantN << endl;
	}

	if (leaves>0){
		auxL = new ulong[cantN];	// these are auxiliary vectors to facility the compute of intervals in each internal nodes. These will be delete later.
		auxR = new ulong[cantN];

		// Step 1: process the n bits in groups of size S, and create the leaves in each array:	intervals of relatives excess per block, in MinLeaf and MaxLeaf,
		//         and the relative sum of excess per block, in ELeaf
		// 'i' is for bits in P, 'cont' is for count the size of each block, and 'node' is the current leaf position to set min-max values...
		//cout << "Step 1: process the n bits in groups of size s, and create the leaves in each array..." << endl;
		j = firstLeaf; // this is for auxiliary vectors
		ulong bitBottom = leavesBottom*Srmq;
		for (i=node=cont=0; i<nBin; i++){
			if (i == bitBottom){
				node = leavesBottom;	// we move up one level
				j = cantIN;
			}
			if(cont==0){
				auxL[j] = i;
				position = j;
				// set left boundaries, when the index j is the first child.
				while(position && (position+1)%2==0){
					father = (position+1)>>1;
					if ((position+1)%2 > 1)
						father++;
					father--;
					auxL[father] = i;
					position = father;
				}
				cont=1;
			}else{
				cont++;
				if(cont == Srmq){
					auxR[j] = i;
					position = j;
					// set right boundaries, when the index j is the last child. The last child always has rest == 1.
					while(position && (position+1)%2==1){
						father = (position+1)>>1;
						if ((position+1)%2 > 1)
							father++;
						father--;
						auxR[father] = i;
						position = father;
					}
					node++;
					j++;
					cont = 0;
				}
			}
		}
		if (cont)
			auxR[node] = i;
	}

	// Step 2: create arrays of relative and super blocks...
	//cout << "Step 2: create arrays of relative and super blocks..." << endl;
	createTables();
	if (TRACE && leaves>0){
		cout << endl << "min-max Intervals..." << endl;
		for (i=0; i<cantN; i++)
			cout << "[" << auxL[i] << "," << auxR[i] << "] ";
		cout << endl;
	}

	MIN_BCK = 0;
	if (leaves>0){
		ulong segment;
		long int currSumBck;
		Aux_Bkwd_MinIN = new int[cantIN];
		MAX_B = 0;
		MIN_BCK = 1;
		// Step 3: set the min - max relatives values to each internal node.
		//cout << "Step 3: set the min - max relatives values to each internal node..." << endl;
		for(i=cantIN; i>0; i--){
			child = i<<1;		// child is the second child of node i
			if (child > cantIN){
				if (child > firstLeaf)
					node = child - firstLeaf;
				else
					node = child - cantIN + leavesBottom;
				currSumBck = 0;
				miniBck = 0; // this 'min' correspond to the maximum level in the tree
				for (j=0, rb=N8W64-1; j<N8Srmq; j++){
					segment = (P[((node+1)*Srmq-1-BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
					if (currSumBck + T_MIN_BCK[segment] > miniBck)
						miniBck = currSumBck + T_MIN_BCK[segment];
					currSumBck += T_SUM_BLOCK[segment];
					if (rb == 0) rb=N8W64-1;
					else rb--;
				}
				for (j=0, rb=N8W64-1; j<N8Srmq; j++){
					segment = (P[(node*Srmq-1-BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
					if (currSumBck + T_MIN_BCK[segment] > miniBck)
						miniBck = currSumBck + T_MIN_BCK[segment];
					currSumBck += T_SUM_BLOCK[segment];
					if (rb == 0) rb=N8W64-1;
					else rb--;
				}
			}else{
				miniBck = Aux_Bkwd_MinIN[child];
				currSumBck = sumAtPos(auxR[child]) - sumAtPos(auxL[child]-1);
				if (currSumBck + Aux_Bkwd_MinIN[child-1] > miniBck)
					miniBck = currSumBck + Aux_Bkwd_MinIN[child-1];

			}
			Aux_Bkwd_MinIN[i-1] = miniBck;
			if (miniBck > MIN_BCK)
				MIN_BCK = miniBck;

		}
		delete [] auxL;
		delete [] auxR;
	}

	lgMIN_BCK = ceilingLog64(MIN_BCK+1, 2);
	if (TRACE) cout << "MIN_BCK " << MIN_BCK << ", lgMIN_BCK " << lgMIN_BCK << endl;

	sizeDS = cantIN*lgMIN_BCK/W64;
	if ((cantIN*lgMIN_BCK)%W64)
		sizeDS++;
	Bkwd_MinIN = new ulong[sizeDS];
	sizeDS = sizeDS*sizeof(ulong);
	sizeRMM += sizeDS;
	if (TRACE) cout << " ** size of Bkwd_MinIN[] " << sizeDS << " Bytes" << endl;

	ulong cMIN_BCK = 0; // count
	for(i=0; i<cantIN; i++){
		setNum64(Bkwd_MinIN, cMIN_BCK, lgMIN_BCK, Aux_Bkwd_MinIN[i]);
		cMIN_BCK += lgMIN_BCK;
	}
	if (cantIN)
		delete [] Aux_Bkwd_MinIN;

	if (TRACE) cout << " ** Total RMQRMM64 size: " << sizeRMM << " Bytes = " << (float)sizeRMM/(1024.0*1024.0) << " MB." << endl;
}

void RMQRMM64::createTables(){
	ulong sizeDS;
	ulong i, j, jS, jR, cont, rb, segment, posMin;
	long int sum, sumBlock;
	int Min;
	lenSB = (nP/Srmq)/SuBrmq + 1;
	//bitsSuB = Srmq*SuBrmq;
	ulong *AuxTSBlock = new ulong[lenSB];

	TRBlock = new char[lenSB];
	sizeDS = lenSB*sizeof(char);
	sizeRMM += sizeDS;
	if (TRACE) cout << " ** size of TRBlock[] " << sizeDS << " Bytes" << endl;

	ulong lenWRel = lenSB/W64;
	if (lenSB%W64)
		lenWRel++;
	Bfull = new ulong[lenWRel];
	sizeDS = lenWRel*sizeof(ulong);
	sizeRMM += sizeDS;
	if (TRACE) cout << " ** size of Bfull[] " << sizeDS << " Bytes" << endl;

	j = leaves/W64;
	if (j%W64)
		j++;
	TPMinB = new uchar[leaves];		// positions explicitly
	sizeDS = leaves*sizeof(uchar);
	sizeRMM += sizeDS;
	if (TRACE) cout << " ** size of TPMinB[] " << sizeDS << " Bytes" << endl;

	int *AuxTMinB = new int[leaves];
	jR = cont = sum = 0;
	jS = 1;
	AuxTSBlock[0] = TRBlock[0] = MAX_SupB = 0;
	// here.... always sum >= 0, because is a BP sequence
	for (i=0; jS<lenSB; i++){
		rb = sumBlock = 0;
		for (j=0; j<N8Srmq; j++){
			segment = (P[(i*Srmq+BSrmq*j)/W64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			//printBitsUlong(segment);cout<<endl;
			sumBlock += T_SUM_BLOCK[segment];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}

		sum += sumBlock;
		if(cont==0){
			if (sumBlock==-Srmq || sumBlock==Srmq){
				setBit64(Bfull,jR);
				if(sumBlock==-Srmq)
					TRBlock[jR] = 0;
				else
					TRBlock[jR] = 1;
			}else{
				cleanBit64(Bfull,jR);
				TRBlock[jR] = (char)(sumBlock>>1);
			}
			jR++;
		}
		cont++;
		if (cont == SuBrmq){
			AuxTSBlock[jS] = sum;
			if (sum > (int)MAX_SupB)
				MAX_SupB = sum;
			jS++;
			cont = 0;
		}
	}

	if((Srmq+((lenSB-1)<<bitsSuB)-1)<nP){
		rb = sumBlock = 0;
		for (j=0; j<N8Srmq; j++){
			segment = (P[((jR<<bitsSuB)+BSrmq*j)/W64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			sumBlock += T_SUM_BLOCK[segment];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}

		if (sumBlock==-Srmq || sumBlock==Srmq){
			setBit64(Bfull,jR);
			if(sumBlock==-Srmq)
				TRBlock[jR] = 0;
			else
				TRBlock[jR] = 1;
		}else{
			cleanBit64(Bfull,jR);
			TRBlock[jR] = (char)(sumBlock>>1);
		}
	}

	if (lenSB == 1){
		rank1_Bin = 0;
		for (j=0; j<nBin; j++){
			if(readBit64(P, j))
				rank1_Bin++;
		}
	}else
		rank1_Bin = ((nBin-sum)>>1) + sum;
	sum = 0;
	for (i=0; i<h; i++){
		cont = (ulong)pow(2, (double)(i));
		sum += cont;
	}

	if (leaves){
		int maxMin = -1*(pow(2,31)-1);
		for (i=1; i<=leaves; i++){
			sumBlock = posMin = 0;
			Min = -1*(pow(2,31)-1);

			cont = i*Srmq-1;
			for (j=1; j<=Srmq; j++, cont--){
				if (readBit64(P, cont)){
					sumBlock++;
					if(Min < sumBlock){
						Min = sumBlock;
						posMin = Srmq-cont%Srmq-1;
					}
				}else sumBlock--;
			}
			TPMinB[i-1] = posMin;
			AuxTMinB[i-1] = Min;
			if (maxMin < Min)
				maxMin = Min;
		}
		if (maxMin < 0){
			cout << "ERR... maxMin = " << maxMin << " < 0" << endl;
			exit(0);
		}
		MAX_B = maxMin;
		lgMAX_SupB = ceilingLog64(MAX_SupB+1, 2);		// include one bit for sign
	}else{
		MAX_B = 0;
		lgMAX_SupB = 1;
	}

	if(MAX_B == 0)
		lgMAX_B = 1;
	else
		lgMAX_B = ceilingLog64(MAX_B+1, 2);

	cont = leaves*lgMAX_B/W64;
	if ((leaves*lgMAX_B)%W64)
		cont++;
	TMinB = new ulong[cont];
	sizeDS = cont*sizeof(ulong);
	sizeRMM += sizeDS;
	if (TRACE) cout << " ** size of TMinB[] " <<  sizeDS << " Bytes" << endl;

	for (i=cont=0; i<leaves; i++){
		if (AuxTMinB[i] < 0)
			setNum64(TMinB, cont, lgMAX_B, 0);
		else
			setNum64(TMinB, cont, lgMAX_B, (ulong)AuxTMinB[i]);
		cont += lgMAX_B;
	}
	if (leaves)
		delete []AuxTMinB;

	cont = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		cont++;
	TSBlock = new ulong[cont];
	sizeDS = cont*sizeof(ulong);
	sizeRMM += sizeDS;
	if (TRACE) cout << " ** size of TSBlock[] " <<  sizeDS << " Bytes" << endl;

	for (i=cont=0; i<lenSB; i++, cont+=lgMAX_SupB)
		setNum64(TSBlock, cont, lgMAX_SupB, AuxTSBlock[i]);
	if (lenSB)
		delete [] AuxTSBlock;

	if (TRACE){
		cout << "MAX_SupB " << MAX_SupB << ", lgMAX_SupB " << lgMAX_SupB << ", MAX_B " << MAX_B << ", lgMAX_B " << lgMAX_B << endl;
		cout << "TSBlock[1.." <<lenSB<< "]..." << endl;
		for (i=0; i<lenSB; i++)
			cout << getNum64(TSBlock, i*lgMAX_SupB, lgMAX_SupB) << " ";
		cout << endl;
		cout << "TRBlock[1.." <<lenSB<< "]..." << endl;
		for (i=0; i<lenSB; i++)
			cout << (int)TRBlock[i] << " ";
		cout << endl;
		cout << "Bfull[1.." <<lenSB<< "]..." << endl;
		for (i=0; i<lenSB; i++)
			cout << readBit64(Bfull, i);
		cout << endl;
		cout << "TPMinB[1.." <<leaves<< "]..." << endl;
		for (i=0; i<leaves; i++)
			cout << (uint)TPMinB[i] << " ";
		cout << endl;
		cout << "TMinB[1.." <<leaves<< "]..." << endl;
		for (i=0; i<leaves; i++)
			cout << getNum64(TMinB, i*lgMAX_B, lgMAX_B) << " ";
		cout << endl;
	}
}

// *******************************************************************************
// ********************************* BASIC OPERATIONS ****************************

// for 0 <= i < n
ulong RMQRMM64::binRank_1(ulong i){
	return (sumAtPos(i)+i+1)>>1;
}

ulong RMQRMM64::binSelect_1(ulong i){
	ulong blk, blk2, rank_B, raAux, q, rb, curr;

	// search on super blocks...
	blk = (i>>bitsSuB)<<1;  // proportional search considering that it is a sequence balanced of 1's and 0's
	if (blk){
		raAux = (blk<<bitsSuB + getNum64(TSBlock, blk*lgMAX_SupB, lgMAX_SupB))>>1;
		if(raAux<i){
			rank_B = raAux;
			blk2 = blk+1;
			raAux = (blk2<<bitsSuB + getNum64(TSBlock, blk2*lgMAX_SupB, lgMAX_SupB))>>1;
			while(raAux<i){
				rank_B = raAux;
				blk = blk2;
				blk2++;
				raAux = (blk2<<bitsSuB + getNum64(TSBlock, blk2*lgMAX_SupB, lgMAX_SupB))>>1;
			}
		}else{
			while(raAux >= i){		// I exceeded! come back...
				blk--;
				raAux = (blk<<bitsSuB + getNum64(TSBlock, blk*lgMAX_SupB, lgMAX_SupB))>>1;
			}
			rank_B = raAux;
		}

		// search in blocks...
		if(readBit64(Bfull, blk)){
			if(TRBlock[blk])
				raAux = rank_B + Srmq;
			else
				raAux = rank_B;
		}else
			raAux = rank_B + TRBlock[blk] + SrmqM;

		curr = SrmqD*blk;
		if(raAux<i){
			rank_B = raAux;
			curr += Srmq;
		}
	}else{
		if(readBit64(Bfull, 0)){
			if(TRBlock[0])
				raAux = Srmq;
			else
				raAux = 0;
		}else
			raAux = SrmqM + TRBlock[0];

		if(raAux<i){
			rank_B = raAux;
			curr=Srmq;
		}else
			rank_B = curr = 0;
	}

	// search in the leaves... go from the left... it is possible go from the right !!
	q = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
	raAux = rank_B + ((BSrmq + T_SUM_BLOCK[q])>>1);
	rb = 0;
	while (raAux < i){
		rank_B = raAux;
		curr+=BSrmq;
		rb++;
		if (rb == N8W64)
			rb=0;
		q = (P[curr>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
		raAux += (BSrmq + T_SUM_BLOCK[q])>>1;
	}
	for (; rank_B<i; curr++){
		if (readBit64(P, curr))
			rank_B++;
	}

	return curr-1;
}

// for 0 <= i < n
ulong RMQRMM64::rank_1(ulong i){
	if(i < nBin-1)
		return binRank_1(i);
	else{
		if (nBin-1 == i) return rank1_Bin;
		if(i >= nP-1) return nP>>1;

		ulong blk=i>>bitsSuB;
		ulong x = blk<<bitsSuB;
		ulong rank = (x + getNum64(TSBlock, blk*lgMAX_SupB, lgMAX_SupB))>>1;
		if((i-x) > Srmq){
			if(readBit64(Bfull, blk)){
				if(TRBlock[blk])
					rank += Srmq;
			}else
				rank += SrmqM + TRBlock[blk];
			x+=Srmq;
		}
		ulong b, rest, q;
		while(x+BrmqMOne <= i){
			b = x>>BW64;
			rest = (x+BSrmq)%W64;
			q = (P[b] >> (W64-rest)) & 0xff;
			rank += __popcount_tab[q];
			x += BSrmq;
		}

		// check last segment (< S) bit by bit...
		while (x<=i){
			if(readBit64(P,x))
				rank++;
			x++;
		}

		return rank;
	}
}

ulong RMQRMM64::select_1_new(ulong i){
	ulong nxt, curr=0, pos, l=0, r, m, ones=0;

	if (lenSB>1){
		r=lenSB;
		while (l<=r){
			m = l+((r-l)>>1);
			nxt = ((m<<bitsSuB)+getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
			while (nxt<i){
				ones = nxt;
				curr = m;
				if (l<r){
					l = m+1;
					m = l+((r-l)>>1);
					nxt = ((m<<bitsSuB)+getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
				}else
					break;
			}
			r = m-1;
			if (l<=r){
				m = l+((r-l)>>1);
				nxt = ((m<<bitsSuB)+getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
				while (nxt>i && l<m){
					r = m-1;
					m = l+((r-l)>>1);
					nxt = ((m<<bitsSuB)+getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
				}
				if(nxt<i){
					ones = nxt;
					curr = m;
				}
				l = m+1;
			}else
				break;
		}


		if(readBit64(Bfull, curr)){
			if(TRBlock[m])
				nxt = ones+Srmq;
			else
				nxt = ones;
		}else
			nxt = ones+SrmqM+TRBlock[curr];
		curr = (curr<<bitsSuB)+BrmqMOne;
		if (nxt < i){
			ones = nxt;
			curr += Srmq;
		}
	}else
		curr = BrmqMOne;

	r = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
	nxt = ones+__popcount_tab[r];
	for(l=1; nxt<i; l++){
		if (l==N8W64)
			l=0;
		curr+=BSrmq;
		r = (P[curr>>BW64] & RMMMasks[l]) >> (W64m8-BSrmq*l);
		nxt += __popcount_tab[r];
	}
	ones = nxt-__popcount_tab[r];
	pos=curr-BrmqMOne;
	for (; ones<i; pos++){
		if (readBit64(P, pos))
			ones++;
	}

	return pos-1;
}

ulong RMQRMM64::select_1(ulong i){
	ulong nxt, curr=0, l=0, r, m, ones=0;

	if (i <= rank1_Bin && lenSB>1){
		r=lenSB;
		while (l<=r){
			m = l+((r-l)>>1);
			nxt = ((m<<bitsSuB)+getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
			while (nxt<i){
				ones = nxt;
				curr = m;
				if (l<r){
					l = m+1;
					m = l+((r-l)>>1);
					nxt = ((m<<bitsSuB)+getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
				}else
					break;
			}
			r = m-1;
			if (l<=r){
				m = l+((r-l)>>1);
				nxt = ((m<<bitsSuB)+getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
				while (nxt>i && l<m){
					r = m-1;
					m = l+((r-l)>>1);
					nxt = ((m<<bitsSuB)+getNum64(TSBlock, m*lgMAX_SupB, lgMAX_SupB))>>1;
				}
				if(nxt<i){
					ones = nxt;
					curr = m;
				}
				l = m+1;
			}else
				break;
		}

		if(readBit64(Bfull, curr)){
			if(TRBlock[m])
				nxt = ones+Srmq;
			else
				nxt = ones;
		}else
			nxt = ones+SrmqM+TRBlock[curr];
		curr = (curr<<bitsSuB)+BrmqMOne;
		if (nxt < i){
			ones = nxt;
			curr += Srmq;
		}
	}else{
		curr = nBin+BrmqMOne;
		ones = rank1_Bin;
	}
	r = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
	ones+=__popcount_tab[r];
	for(l=1; ones<i; l++){
		if (l==N8W64)
			l=0;
		curr+=BSrmq;
		r = (P[curr>>BW64] & RMMMasks[l]) >> (W64m8-BSrmq*l);
		ones += __popcount_tab[r];
	}
	ones-=__popcount_tab[r];
	curr-=BrmqMOne;
	for (; ones<i; curr++){
		if (readBit64(P, curr))
			ones++;
	}

	return curr-1;
}

ulong RMQRMM64::select_1_old(ulong i){
	ulong pos = 0;

	if(i <= rank1_Bin)
		return binSelect_1(i);
	else{
		ulong rb, l, q, curr;
		ulong rank_A, rank_B = rank1_Bin;

		curr = nBin+BrmqMOne;
		q = (P[curr>>BW64] & RMMMasks[0]) >> W64m8;
		rank_A = rank_B + ((BSrmq + T_SUM_BLOCK[q])>>1);
		rb = l = 0;
		while (rank_A < i){
			rank_B = rank_A;
			rb++;
			if (rb == N8W64)
				rb=0;
			curr+=BSrmq;
			q = (P[curr>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			rank_A += (BSrmq + T_SUM_BLOCK[q])>>1;
		}
		pos=curr-BrmqMOne;
		for (; rank_B<i; pos++){
			if (readBit64(P, pos))
				rank_B++;
		}
	}

	return pos-1;
}

// give the excess from 0 to pos
long int RMQRMM64::sumAtPos(long int pos){
	ulong rb, q, i, l, k, j=pos+1;
	ulong blk = j>>PotSrmq;
	ulong sup = blk/SuBrmq;
	long int sum = getNum64(TSBlock, sup*lgMAX_SupB, lgMAX_SupB);

	if(blk%2){
		if(readBit64(Bfull,sup)){
			if(TRBlock[sup])
				sum += Srmq;
			else
				sum -= Srmq;
		}else
			sum += TRBlock[sup] << 1;
	}

	if (j>=BSrmq){
		for (rb=0, l=blk<<PotSrmq, k=j-BSrmq; l<=k; l+=BSrmq){
			q = (P[l>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			sum += T_SUM_BLOCK[q];
			if (rb == N8W64-1) rb=0;
			else rb++;
		}
	}
	if (j%BSrmq){
		for (i=pos-j%BSrmq+1; i<j; i++){
			if (readBit64(P, i))
				sum++;
			else
				sum--;
		}
	}
	return sum;
}

// *******************************************************************************
// ****************************** TREE OPERATIONS ********************************

// give the excess of the internal node 'node=preorder+1' that has a distance 'dist' to the tree's depth
long int RMQRMM64::computeSumOfNode(ulong node, ulong dist){
	long int sum = 0;
	ulong ini=node<<dist, end=(node+1)<<dist;

	if(ini>cantN){
		ini>>=1;
		end>>=1;
	}
	if (ini >= firstLeaf)
		ini -= firstLeaf+1;
	else
		ini = ini-cantIN+leavesBottom-1;

	if (end >= cantN){
		end>>=1;

		if(end < firstLeaf)
			end = leavesBottom+end-cantIN-1;
		else
			end = leaves;
	}else{
		if (end > firstLeaf)
			end -= firstLeaf+1;
		else
			end = end-cantIN+leavesBottom-1;
	}
	if(ini > end) end = leaves;

	sum += sumAtPos((end<<PotSrmq)-1) - sumAtPos((ini<<PotSrmq)-1);
	return sum;
}

// give the number of leaves of this node 'node=preorder+1' that has a distance 'dist' to the tree's depth
ulong RMQRMM64::computeLeavesOfNode(ulong node, ulong dist){
	ulong ini=node<<dist, end=(node+1)<<dist;

	if(ini>cantN){
		ini>>=1;
		end>>=1;
	}
	if (ini >= firstLeaf)
		ini -= firstLeaf+1;
	else
		ini = ini-cantIN+leavesBottom-1;

	if (end >= cantN){
		end>>=1;

		if(end < firstLeaf)
			end = leavesBottom+end-cantIN-1;
		else
			end = leaves;
	}else{
		if (end > firstLeaf)
			end -= firstLeaf+1;
		else
			end = end-cantIN+leavesBottom-1;
	}
	if(ini > end) end = leaves;

	return end-ini;
}

// search the minimum in a block (sequential search)
void RMQRMM64::search_min_block(ulong x1, ulong x2, long int *min, long int *curSum, ulong *position){
	int Min, sum, auxMin;
	uint rest, q;
	ulong b, len = x2-x1+1, posMin = *position;

	Min = *min;
	sum = *curSum;
	while(len > BrmqMOne){
		b = x2>>BW64;
		rest = (x2+1)%W64;
		if (b == (x2-BrmqMOne)>>BW64)				// x and (x-7) are in the same word...
			q = (P[b] >> (W64-rest)) & 0xff;
		else
			q = ((P[b-1] << rest) | (P[b] >> (W64-rest))) & 0xff;
		auxMin = T_MIN_BCK[q];
		if (Min < sum+auxMin){
			Min = sum+auxMin;
			posMin = x2-T_BCK_D[q][auxMin-1];
		}
		sum += T_SUM_BLOCK[q];
		len -= BSrmq;
		x2 -= BSrmq;
	}
	if (len){
		// check last segment (len < 8) bit by bit...
		while (len>=1){
			if(readBit64(P,x2)){
				sum++;
				if (sum > Min){
					Min = sum;
					posMin = x2;
				}
			}else sum--;
			x2--;
			len--;
		}
	}
	*curSum = sum;
	*min = Min;
	*position = posMin;
}

ulong RMQRMM64::rmqi(ulong i, ulong j){
	long int min, sum;
	min = sum = 0;
	if(j < nBin)
		return rmqi_rmm(i, j, &min, &sum, j);
	else{
		ulong posMin = j;

		if(i >= nBin){
			search_min_block(i, j, &min, &sum, &posMin);
			return posMin;
		}else{
			if (nBin){
				search_min_block(nBin, j, &min, &sum, &posMin);
				return rmqi_rmm(i, nBin-1, &min, &sum, posMin);
			}
		}
	}
	return 0;
}

// return the position of the open parenthesis closet to the root between i and j
ulong RMQRMM64::rmqi_rmm(ulong x1, ulong x2, long int *min, long int *currSum, ulong posMin){
	ulong lLeaf, aLeaf, node, nodeMin;
	ulong height, group, distNodeMin, dist = 1;
	long int mini, sumMin;
	bool isInNode = false;
	bool isMinInLeaf = false;

	lLeaf = x1>>PotSrmq;
	node = x2>>PotSrmq;

	// [1]- search the minimum inside the rightmost block...
	if (lLeaf == node){
		// here is the answer
		search_min_block(x1, x2, min, currSum, &posMin);
		return posMin;
	}else
		search_min_block((node<<PotSrmq), x2, min, currSum, &posMin);
	long int Min = *min, sum = *currSum;

	// [2]- We looking for in the next leaf
	if (node%2){	// this is a right leaf
		node--;
		if (node == lLeaf){		// here is the answer
			search_min_block(x1, ((node+1)<<PotSrmq)-1, min, currSum, &posMin);
			return posMin;
		}
		mini = getNum64(TMinB, node*lgMAX_B, lgMAX_B);
		if (Min < sum + mini){
			Min = sum + mini;
			posMin = ((node+1)<<PotSrmq) - TPMinB[node] -1;
			isInNode = isMinInLeaf = true;
		}
		if(readBit64(Bfull, (node>>1))){
			if(TRBlock[node>>1])
				sum += Srmq;
			else
				sum -= Srmq;
		}else
			sum += TRBlock[node>>1]<<1;
	}
	if (node == 0){
		*min = Min;
		return posMin;
	}
	aLeaf = node;

	// [3]- We climb recomputing the Min until the position i.
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
	node>>= 1;
	node--;
	dist++;
	group = computeLeavesOfNode(node+1, dist);
	while (lLeaf < aLeaf-group){                  /////// group > aLeaf && ...
		aLeaf -= group;
		mini = sum+getNum64(Bkwd_MinIN, node*lgMIN_BCK, lgMIN_BCK);
		if (Min < mini){
			Min = mini;
			isInNode = true;
			isMinInLeaf = false;
			nodeMin = node;
			sumMin = sum;
			distNodeMin = dist;
		}
		sum += computeSumOfNode(node+1, dist);
		if (node%2==0)
			node--;
		else{
			node>>=1;		// go to my uncle
			node--;
			height--;
			dist++;
			group<<=1;
		}
		group = computeLeavesOfNode(node+1, dist);
	}

	// [4]- We move down recomputing the Min and reach to the left boundaries leaf.
	node++;
	node<<= 1;
	dist--;
	while (node < cantIN){
		group = computeLeavesOfNode(node+1, dist);
		if(lLeaf < aLeaf-group){
			aLeaf -= group;
			mini = sum+getNum64(Bkwd_MinIN, node*lgMIN_BCK, lgMIN_BCK);
			if (Min < mini){
				Min = mini;
				isInNode = true;
				isMinInLeaf = false;
				nodeMin = node;
				sumMin = sum;
				distNodeMin = dist;
			}
			sum += computeSumOfNode(node+1, dist);
			node--;
		}else{
			node++;
			node<<= 1;
			dist--;
		}
	}
	if (node > firstLeaf)			// at this point, 'node' is a leaf; node >= firstLeaf
		node -= firstLeaf;
	else
		node += leavesBottom - cantIN;

	// [5] looking for in the leaves of the leftmost side
	if (node == lLeaf){
		if (isInNode){
			mini = Min;
			search_min_block(x1, ((node+1)<<PotSrmq)-1, &Min, &sum, &posMin);
			if (mini < Min)
				return posMin;
		}else{
			search_min_block(x1, ((node+1)<<PotSrmq)-1, &Min, &sum, &posMin);
			return posMin;
		}
	}else{
		if (isInNode){
			aLeaf = node>>1;
			mini = getNum64(TMinB, node*lgMAX_B, lgMAX_B);
			if (Min < sum + mini){
				Min = sum + mini;
				posMin = ((node+1)<<PotSrmq) - TPMinB[node] -1;
				sum += getNum64(TSBlock, (aLeaf+1)*lgMAX_SupB, lgMAX_SupB) - getNum64(TSBlock, aLeaf*lgMAX_SupB, lgMAX_SupB);
				if(readBit64(Bfull, aLeaf)){
					if(TRBlock[aLeaf])
						sum -= Srmq;
					else
						sum += Srmq;
				}else
					sum -= TRBlock[aLeaf]<<1;
				search_min_block(x1, (node<<PotSrmq)-1, &Min, &sum, &posMin);
				return posMin;
			}

			sum += getNum64(TSBlock, (aLeaf+1)*lgMAX_SupB, lgMAX_SupB) - getNum64(TSBlock, aLeaf*lgMAX_SupB, lgMAX_SupB);
			if(readBit64(Bfull, aLeaf)){
				if(TRBlock[aLeaf])
					sum -= Srmq;
				else
					sum += Srmq;
			}else
				sum -= TRBlock[aLeaf]<<1;

			mini = Min;
			search_min_block(x1, (node<<PotSrmq)-1, &mini, &sum, &posMin);
			if (mini > Min)
				return posMin;
		}else{
			aLeaf = node>>1;
			mini = getNum64(TMinB, node*lgMAX_B, lgMAX_B);
			if (Min < sum + mini){
				Min = sum + mini;
				posMin = ((node+1)<<PotSrmq) - TPMinB[node] -1;
			}
			sum += getNum64(TSBlock, (aLeaf+1)*lgMAX_SupB, lgMAX_SupB) - getNum64(TSBlock, aLeaf*lgMAX_SupB, lgMAX_SupB);
			if(readBit64(Bfull, aLeaf)){
				if(TRBlock[aLeaf])
					sum -= Srmq;
				else
					sum += Srmq;
			}else
				sum -= TRBlock[aLeaf]<<1;
			search_min_block(x1, (node<<PotSrmq)-1, &Min, &sum, &posMin);
			return posMin;
		}
	}

	if(isMinInLeaf)
		return posMin;
	// [6] Here the minimum value is in a block which descent of the internal node 'nodeMin',
	// then we must descent in order to find its position.
	nodeMin++;
	nodeMin<<= 1;
	distNodeMin--;
	sum = sumMin;
	while (nodeMin < cantIN){
		mini = getNum64(Bkwd_MinIN, nodeMin*lgMIN_BCK, lgMIN_BCK);
		sumMin = computeSumOfNode(nodeMin+1, distNodeMin);
		if (Min == sum+mini){   // is it ? then go to the right without updating the Min
			nodeMin++;
			nodeMin<<= 1;
			distNodeMin--;
		}else{					// else, go to the left updating the sum
			nodeMin--;
			sum += sumMin;
		}
	}
	if (nodeMin > firstLeaf)
		nodeMin -= firstLeaf;
	else
		nodeMin += leavesBottom - cantIN;

	// Here, I get the right leaf...
	aLeaf = nodeMin>>1;
	sumMin = getNum64(TSBlock, (aLeaf+1)*lgMAX_SupB, lgMAX_SupB) - getNum64(TSBlock, aLeaf*lgMAX_SupB, lgMAX_SupB);
	if(readBit64(Bfull, aLeaf)){
		if(TRBlock[aLeaf-1])
			sumMin -= Srmq;
		else
			sumMin += Srmq;
	}else
		sumMin -= TRBlock[aLeaf]<<1;

	if ((long int)getNum64(TMinB, nodeMin*lgMAX_B, lgMAX_B) >= (long int)getNum64(TMinB, (nodeMin-1)*lgMAX_B, lgMAX_B)+sumMin)
		return ((nodeMin+1)<<PotSrmq)-TPMinB[nodeMin]-1;
	else
		return (nodeMin<<PotSrmq)-TPMinB[nodeMin-1]-1;
}

ulong RMQRMM64::queryRMQ(ulong i, ulong j){
	return rank_1(rmqi(select_1(i+2),select_1(j+2)))-2;
}

uint RMQRMM64::getSize(){
	return sizeRMM;
}

void RMQRMM64::saveDS(char *fileName){
	cout << "Save data structure in " << fileName << endl;
	ofstream os (fileName, ios::binary);
	cout << "   Data structure size: " << sizeRMM << endl;

	if(TRACE){
		cout << "Variables load: " << endl;
		cout << "nP " << nP << endl;
		cout << "nW " << nW << endl;
		cout << "rank1_Bin " << rank1_Bin << endl;
		cout << "nBin " << nBin << endl;
		cout << "cantN " << cantN << endl;
		cout << "cantIN " << cantIN << endl;
		cout << "leaves " << leaves << endl;
		cout << "leavesBottom " << leavesBottom << endl;
		cout << "firstLeaf " << firstLeaf << endl;
		cout << "lenSB " << lenSB << endl;
		cout << "lenLB " << lenLB << endl;
		cout << "h " << h << endl;
		cout << "MAX_BCK " << MAX_BCK << endl;
		cout << "MAX_B " << MAX_B << endl;
		cout << "lgMAX_B " << lgMAX_B << endl;
		cout << "MAX_SupB " << MAX_SupB << endl;
		cout << "lgMAX_SupB " << lgMAX_SupB << endl;
		cout << "lgMIN_BCK " << lgMIN_BCK << endl;
		cout << "MIN_BCK " << MIN_BCK << endl;
	}

	os.write((const char*)&nP, sizeof(ulong));
	os.write((const char*)&nW, sizeof(ulong));
	os.write((const char*)&rank1_Bin, sizeof(ulong));
	os.write((const char*)&nBin, sizeof(ulong));
	os.write((const char*)&cantN, sizeof(ulong));
	os.write((const char*)&cantIN, sizeof(ulong));
	os.write((const char*)&leaves, sizeof(ulong));
	os.write((const char*)&leavesBottom, sizeof(ulong));
	os.write((const char*)&firstLeaf, sizeof(ulong));
	os.write((const char*)&lenSB, sizeof(ulong));
	os.write((const char*)&lenLB, sizeof(uint));
	os.write((const char*)&h, sizeof(uint));
	os.write((const char*)&MAX_BCK, sizeof(uint));
	os.write((const char*)&MAX_B, sizeof(uint));
	os.write((const char*)&lgMAX_B, sizeof(uint));
	os.write((const char*)&MAX_SupB, sizeof(uint));
	os.write((const char*)&lgMAX_SupB, sizeof(uint));
	os.write((const char*)&lgMIN_BCK, sizeof(uint));
	os.write((const char*)&MIN_BCK, sizeof(int));

	ulong sizeDT = 10*sizeof(ulong) + 9*sizeof(uint) + sizeof(int);
	sizeDT +=  2*512 + 2048;											// size for T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[]
	//if(TRACE)
		cout << " .- T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[] + Variables " << sizeDT << " Bytes" << endl;

	ulong size = nP >> BW64;
	if (nP % W64)
		size++;
	os.write((const char*)P, size*sizeof(ulong));				// save P[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- P[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = cantIN*lgMIN_BCK/W64;
	if ((cantIN*lgMIN_BCK)%W64)
		size++;
	os.write((const char*)Bkwd_MinIN, size*sizeof(ulong));		// save Bkwd_MinIN[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- Bkwd_MinIN[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		size++;
	os.write((const char*)TSBlock, size*sizeof(ulong));			// save TSBlock[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- TSBlock[] " << size*sizeof(ulong) << " Bytes" << endl;

	os.write((const char*)TRBlock, lenSB*sizeof(char));			// save TRBlock[]
	sizeDT += lenSB*sizeof(char);
	if(TRACE) cout << " .- TRBlock[] " << lenSB*sizeof(char) << " Bytes" << endl;

	size = lenSB/W64;
	if (lenSB%W64)
		size++;
	os.write((const char*)Bfull, size*sizeof(ulong));			// save Bfull[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- Bfull[] " << size*sizeof(ulong) << " Bytes" << endl;

	os.write((const char*)TPMinB, leaves*sizeof(uchar));		// save TPMinB[]  = new uchar[leaves]
	sizeDT += leaves*sizeof(uchar);
	if(TRACE) cout << " .- TPMinB[] " << leaves*sizeof(uchar) << " Bytes" << endl;

	size = leaves*lgMAX_B/W64;
	if ((leaves*lgMAX_B)%W64)
		size++;
	os.write((const char*)TMinB, size*sizeof(ulong));							// save TMinB[]
	sizeDT += size*sizeof(ulong);
	if(TRACE) cout << " .- TMinB[] " << size*sizeof(ulong) << " Bytes" << endl;

	os.close();
	cout << "   Total bytes saved from data structure: " << sizeDT << endl;
}

void RMQRMM64::loadDS(char *fileName){
	cout << " Load data structure from " << fileName << endl;
	ifstream is(fileName, ios::binary);

	is.read((char*)&nP, sizeof(ulong));
	is.read((char*)&nW, sizeof(ulong));
	is.read((char*)&rank1_Bin, sizeof(ulong));
	is.read((char*)&nBin, sizeof(ulong));
	is.read((char*)&cantN, sizeof(ulong));
	is.read((char*)&cantIN, sizeof(ulong));
	is.read((char*)&leaves, sizeof(ulong));
	is.read((char*)&leavesBottom, sizeof(ulong));
	is.read((char*)&firstLeaf, sizeof(ulong));
	is.read((char*)&lenSB, sizeof(ulong));
	is.read((char*)&lenLB, sizeof(uint));
	is.read((char*)&h, sizeof(uint));
	is.read((char*)&MAX_BCK, sizeof(uint));
	is.read((char*)&MAX_B, sizeof(uint));
	is.read((char*)&lgMAX_B, sizeof(uint));
	is.read((char*)&MAX_SupB, sizeof(uint));
	is.read((char*)&lgMAX_SupB, sizeof(uint));
	is.read((char*)&lgMIN_BCK, sizeof(uint));
	is.read((char*)&MIN_BCK, sizeof(int));

	if(TRACE){
		cout << "Variables load: " << endl;
		cout << "nP " << nP << endl;
		cout << "nW " << nW << endl;
		cout << "rank1_Bin " << rank1_Bin << endl;
		cout << "nBin " << nBin << endl;
		cout << "cantN " << cantN << endl;
		cout << "cantIN " << cantIN << endl;
		cout << "leavesBottom " << leavesBottom << endl;
		cout << "leaves " << leaves << endl;
		cout << "firstLeaf " << firstLeaf << endl;
		cout << "lenSB " << lenSB << endl;
		cout << "lenLB " << lenLB << endl;
		cout << "h " << h << endl;
		cout << "MAX_BCK " << MAX_BCK << endl;
		cout << "MAX_B " << MAX_B << endl;
		cout << "lgMAX_B " << lgMAX_B << endl;
		cout << "MAX_SupB " << MAX_SupB << endl;
		cout << "lgMAX_SupB " << lgMAX_SupB << endl;
		cout << "lgMIN_BCK " << lgMIN_BCK << endl;
		cout << "MIN_BCK " << MIN_BCK << endl;
	}

	// size for variables
	sizeRMM = 10*sizeof(ulong) + 9*sizeof(uint) + sizeof(int);
	sizeRMM += 2*512 + 2048;											// size for T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[]
	//if(TRACE)
	cout << " .- T_SUM_BLOCK[] + T_MIN_BCK[] + T_BCK_D[] + Variables " << sizeRMM << " Bytes" << endl;

	ulong sizeDS = nP >> BW64;
	if (nP % W64)
		sizeDS++;
	P = new ulong[sizeDS];
	is.read((char*)P, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(TRACE) cout << " .- P[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	sizeDS = cantIN*lgMIN_BCK/W64;
	if ((cantIN*lgMIN_BCK)%W64)
		sizeDS++;
	Bkwd_MinIN = new ulong[sizeDS];
	is.read((char*)Bkwd_MinIN, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(TRACE) cout << " .- Bkwd_MinIN[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	sizeDS = lenSB*lgMAX_SupB/W64;
	if ((lenSB*lgMAX_SupB)%W64)
		sizeDS++;
	TSBlock = new ulong[sizeDS];
	is.read((char*)TSBlock, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(TRACE) cout << " .- TSBlock[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	TRBlock = new char[lenSB];
	is.read((char*)TRBlock, lenSB*sizeof(char));
	sizeRMM += lenSB*sizeof(char);
	if(TRACE) cout << " .- TRBlock[] " << lenSB*sizeof(char) << " Bytes" << endl;

	sizeDS = lenSB/W64;
	if (lenSB%W64)
		sizeDS++;
	Bfull = new ulong[sizeDS];
	is.read((char*)Bfull, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(TRACE) cout << " .- Bfull[] " << sizeDS*sizeof(ulong) << " Bytes" << endl;

	TPMinB = new uchar[leaves];
	is.read((char*)TPMinB, leaves*sizeof(uchar));
	sizeRMM += leaves*sizeof(uchar);
	if(TRACE) cout << " .- TPMinB[] " << leaves*sizeof(uchar) << " Bytes" << endl;

	sizeDS = leaves*lgMAX_B/W64;
	if ((leaves*lgMAX_B)%W64)
		sizeDS++;
	TMinB = new ulong[sizeDS];
	is.read((char*)TMinB, sizeDS*sizeof(ulong));
	sizeRMM += sizeDS*sizeof(ulong);
	if(TRACE) cout << " .- TMinB[] " << sizeDS*sizeof(sizeDS) << " Bytes" << endl;

	is.close();
	cout << " Data Structure loaded !!" << endl;
}

// *************************** PRINT THE TREE ************************************
void RMQRMM64::printTree(){
	ulong i, j, cant, acum, rb, segment;
	long int currSum, mini;
	uint cMIN_BCK;
	cMIN_BCK = 0;

	cout << "Min-max Binary Tree...[min bck]i(sum)" << endl;
	cant = 1;
	for (acum=j=i=0; i<cantIN; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " | ";
			acum = 0;
		}

		mini = getNum64(Bkwd_MinIN, cMIN_BCK, lgMIN_BCK);
		cout << i << "_[" << mini << "]";
		cMIN_BCK += lgMIN_BCK;
	}

	for (i=leavesBottom; i<leaves; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " | ";
			acum = 0;
		}

		rb = 3;
		currSum = 0;
		mini = 0;
		for (j=0; j<N8Srmq; j++){
			segment = (P[((i+1)*Srmq-1-BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			if (currSum + T_MIN_BCK[segment] > mini)
				mini = currSum + T_MIN_BCK[segment];
			currSum += T_SUM_BLOCK[segment];
			if (rb == 0) rb=N8W64-1;
			else rb--;
		}
		cout << i << "_h[" << mini << "](" << currSum << ") ";
	}
	cout << endl;
	for (i=0; i<leavesBottom; i++, j++, acum++){
		if (j==cant){
			cout << endl;
			j = acum = 0;
			cant *= 2;
		}
		if (acum == 2){
			cout << " | ";
			acum = 0;
		}
		rb = 3;
		currSum = 0;
		mini = 0;
		for (j=0; j<N8Srmq; j++){
			segment = (P[((i+1)*Srmq-1-BSrmq*j)>>BW64] & RMMMasks[rb]) >> (W64m8-BSrmq*rb);
			if (currSum + T_MIN_BCK[segment] > mini)
				mini = currSum + T_MIN_BCK[segment];
			currSum += T_SUM_BLOCK[segment];
			if (rb == 0) rb=N8W64-1;
			else rb--;
		}
		cout << i << "_h[" << mini << "]"<< i << "(" << currSum << ") ";
	}
	cout << endl;
}

RMQRMM64::~RMQRMM64() {

	nP = nW = rank1_Bin = nBin = cantN = cantIN = leaves =
	leavesBottom = firstLeaf = lenSB = lenLB =
	h = MAX_BCK = MAX_B = lgMAX_B = MAX_SupB =
	lgMAX_SupB = lgMIN_BCK = MIN_BCK = sizeRMM = 0;

	//delete [] LastNode;
	delete [] Bkwd_MinIN;
	delete [] TSBlock;
	delete [] TRBlock;
	delete [] Bfull;
	delete [] TPMinB;
	delete [] TMinB;
	cout << " ~ RMQRMM64 destroyed !!" << endl;
}

// *******************************************************************************
// **************************** TEST TREE OPERATIONS *****************************
void RMQRMM64::test_rank_1(){
	long int sumR, rank;
	ulong i,j,k;
	cout << "RMQRMM64::test_rank_1..." << endl;

	/*i=224;
	for (j=sumR=0; j<=i; j++){
		if(readBit64(P, j))
			sumR++;
	}
	rank = rank_1(i);
	cout << "rank_1(" << i <<") = " << rank << ", ? sumR =" << sumR << endl;
	exit(0);*/

	if (nBin>0){
		i = nBin-1;
		for (j=sumR=0; j<=i; j++){
			if(readBit64(P, j))
				sumR++;
		}
		rank = rank_1(i);
		if (sumR != rank){
			cout << "sumR("<<i<<") = " << sumR << endl;
			cout << "rank_1("<<i<<") = " << rank << endl;
			exit(0);
		}
	}

	for (k=0; k<TEST; k++){
		i = (rand() % (nP-1));
		for (j=sumR=0; j<=i; j++){
			if(readBit64(P, j))
				sumR++;
		}
		rank = rank_1(i);
		if (sumR != rank){
			cout << "ERROR !! rank1(" << i << ") = " << rank << " != sumR = " << sumR << endl;
			exit(1);
		}
	}
	cout << "  test_rank_1 OK !!" << endl;
}

void RMQRMM64::test_select_1(){
	ulong i, j, k, sum, pos;

	cout << "RMQRMM64::test_select_1..." << endl;
	/*i=3902;
	for (j=sum=0; sum<i; j++){
		if(readBit64(P, j))
			sum++;
	}
	j--;
	pos = select_1(i);
	cout << "select_1(" << i <<") = " << pos << ", ? j=" << j << endl;
	 */
	//exit(0);

	for (k=0; k<TEST; k++){
		i = (rand() % ((nP>>1)-2)) + 1;
		for (j=sum=0; sum<i; j++){
			if(readBit64(P, j))
				sum++;
		}
		j--;
		pos = select_1(i);
		if (j != pos){
			cout << "ERROR !! select1(" << i << ") = " << pos << " != j = " << j << endl;
			exit(1);
		}
	}
	cout << "  test_select_1 OK !!" << endl;
}

void RMQRMM64::test_sumAtPos(){
	ulong i, j, k;
	long int sum, sumPos;

	cout << "RMQRMM64::test_sumAtPos..." << endl;
	/*i=71;
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
		i = (rand() % (nP-2));
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

void RMQRMM64::test_rmqi(){
	ulong i, j, posMin, posRmqi;
	long int sum, Min, k;

	/*i=0, j=578;
	for (sum=0, Min=0, posMin=i, k=j; k>=(int)i; k--){
		if(readBit64(P, k)){
			sum++;
			if (Min < sum){
				Min = sum;
				posMin = k;
			}
		}else
			sum--;
	}
	cout << "MINIMO(" << i << ", " << j << ") = " << posMin << endl;
	posRmqi = rmqi(i,j);
	if (posRmqi != posMin){
		cout << "ERROR !!... rmqi(" << i << ", " << j << ") = " << posRmqi << " != " << posMin << endl;
		exit(1);
	}*/

	cout << "RMQRMM64::test_rmqi..." << endl;
	i=0;
	ulong sumTest = nP/2;
	while(i<nP){
		while(i<nP && readBit64(P, i)==0)
			i++;

		j=i+1;
		while(j<nP){
			while(j<nP && readBit64(P, j)==0)
				j++;

			if (j<nP){
				for (sum=0, Min=0, posMin=i, k=j; k>=(long int)i; k--){
					if(readBit64(P, k)){
						sum++;
						if (Min < sum){
							Min = sum;
							posMin = k;
						}
					}else
						sum--;
				}

				//if (TRACE) cout << "rmqi(" << i << ", " << j << ") " << endl;
				posRmqi = rmqi(i,j);
				if (posRmqi != posMin){
					cout << "ERROR !!... rmqi(" << i << ", " << j << ") = " << posRmqi << " != " << posMin << endl;
					exit(1);
				}
				//if (TRACE) cout << "= " << posRmqi << endl;
			}
			j+=sumTest;
		}
		i+=sumTest;
	}

	cout << "test_rmqi OK !!" << endl;
}
