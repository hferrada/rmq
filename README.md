rmq library. (DCC'2016 version)
===========

This code creates a data structure to compute, in constant time, the RMQ(i,j) on an array 
of integers. The size for array of considerable length is around of 2.2n bits.

Authors: 
	Hector Ferrada and Gonzalo Navarro
	hferrada@dcc.uchile.cl, gnavarro@dcc.uchile.cl

Description:
	This is an RMQ compressed data structure. The implementation is based on the
	method of Fischer and Heun [1]. The tree representation is a light version of
	Range Min-Max Tree of Sadakane and Navarro [2].
	In order to reduce the size, our rmq uses a simplified version of the Range 
	min-max tree of Navarro and Sadakane [2], we only used the backward minimum 
	array to store the ranges (not maximum and only for forward). You can set 
	"true" the variables "RMQRMM64::TRACE" and "RMQRMM64::RUNTEST" to see the 
	details. Also we included an small example ("test.cpp") to show how to use 
	it and included an small experiment for the time. This works only for 64 bits 
	and supporting long sequences.

Make:
	To make the library just give the command 'make', this will
	create the lib: 'rmqrmm.a'.

Compile:
	To use the library you must compile your program linking 'rmqrmm.a' and include
	the folder "includes/RMQRMM64.h" in your sorce code.
	For example, compiling the file rmqrmmBP.cpp included here:<br />
	g++ rmqrmmBP.cpp -o rmqrmmBP -O3 rmqrmm.a or simply run the command 'make test'. it will create the binary 'rmqrmmBP'. This binary have to recieve 7 parameter:<br />
1.- n: the length of random sequence.<br />
2.- 0/1: to load from a file(0) or create(1) the complete structure.<br />
3.- saveLoadFile: the file (the path will be included) to store the data structure in order that you can load this later.<br />
4.- repetitions: number of repetitions for experiments.<br />
5.- pseudoSorted: 1 indicated that the input array will be encrease pseudo-sorted.<br />
6.- RandomWeight: if the file is in mode pseudo-sorted, then A[i] will be a arondom value in the range [i-weight, i+weight].<br />
7.- resultsFile: the file to store the esperiments's results.<br />

For example, this line execute the code for n=10^4, stores the data in './rmqrmmBP-data.rmq' and the results in rmqrmmBP-Ramdom.txt:<br />
./rmqrmmBP 10000 1 rmqrmmBP-data.rmq 2000000 0 10000 rmqrmmBP-Ramdom.txt


Note about BP-construction:

The current code included in this library computes the BP representation of the isomorphism tree of Fisher and Heum in $O(n)$ time.<br />
The extra space used is linear in words (i.e., $n \log n$ bits).<br />
However, it could easily be reduced to $O(h)$ words as Fisher and Heum porposed, where $h$ is the height of the tree, but the time becomes $O(n^2)$ for all case.<br />

References:
	Please, if you want to include this tool as part of your experiments, in your
	references please you include the reference [3]. 

[1]. J. Fischer and V. Heun. Space-efficient preprocessing schemes for
range minimum queries on static arrays. SIAM Journal on Computing, 40(2):465–492, 2011.

[2]. K. Sadakane and G. Navarro. Fully-Functional Static and Dynamic Succinct Trees. 
ACM Transactions on Algorithms 10(3):article 16, 2014

[3]. H. Ferrada and G. Navarro. Improved Range Minimum Queries.
To appear in Proc. 26th Data Compression Conference (DCC), 2016.
