rmq library
===========

This code creates a data structure to compute the RMQ(i,j) on an array of integers.

Authors: 
	Hector Ferrada and Gonzalo Navarro
	hferrada@dcc.uchile.cl, gnavarro@dcc.uchile.cl

Description:
	This is a RMQ compressed data structure. The implementation is based on the
	method of Fischer and Heun [1]. The tree representation is a light version of
	Range Min-Max Tree of Sadakane and Navarro [2].
	In order to reduce the size, our rmq uses a simplified version of the Range 
	min-max tree of Navarro and Sadakane [2], we only used the backward minimum 
	array to store the ranges (not maximum and only for backward). You can set 
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
	For example, compiling the file test.cpp included here:
		g++ test.cpp -o myRMQ -O3 rmqrmm.a 
		or simply run the command 'make test'. it will create the binary 'myRMQ'.
	This binary have to recieve two parameter: 
		1- n: the length of random sequence 
		2- rmqFile: the file (the path will be included) to store the data structure 
		in order that you can load this later.
	For example, this line execute the code for n=10^6 and store the data in './rmqFile.rmq':
		./myRMQ 1000000 rmqFile.rmq

References:
	Please, if you want to include this tool as part of your experiments, in your
	references include the two papers above. Later, it will appear another 
	publication to replace these ones.

[1]. J. Fischer and V. Heun. Space-efficient preprocessing schemes for range minimum 
queries on static arrays. SIAM '11.

[2]. K. Sadakane and G. Navarro. Fully-functional succinct trees. In Proceedings of 
the Twenty-First Annual ACM-SIAM Symposium on Discrete Algorithms, SODA '10.