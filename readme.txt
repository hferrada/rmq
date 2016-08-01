Authors: 
	Hector Ferrada and Gonzalo Navarro
	hferrada@dcc.uchile.cl, gnavarro@dcc.uchile.cl

Description:
	This is a small library with basic methods usually used in Document Retrieval.
	'DRF' are the initials for Document Retrieval Framewok.
	This include a RMQ compressed data structure. The implementation is based on the
	method of Fischer and Heun [1]. The tree representation is a light version of
	Range Min-Max Tree of Sadakane and Navarro [2].
	In order to reduce the size, our rmq uses a simplified version of the Range 
	min-max tree of Navarro and Sadakane [2], we only used the backward minimum 
	array to store the ranges (not maximum and only for backward). You can set 
	"true" the variables "RMQRMM64::TRACE" and "RMQRMM64::RUNTEST" to see the 
	details. Also we included an small example ("rmqrmmBP.cpp") to show how to use 
	it and included an small experiment for the time. This works only for 64 bits 
	and supporting long sequences.

Make:
	To make the library just give the command 'make', this will
	create the lib: 'rmqrmmBP.a'.

Compile:
	To use the library you must compile your program linking 'rmqrmmBP.a' and include
	the the header "includes/RMQRMM64.h" to buil RMQ compressed structures.
	For example, compiling the file test.cpp included here:
		g++ rmqrmmBP.cpp -o myRMQ -O3 rmqrmmBP.a 
		or simply run the command 'make test'. it will create the binary 'rmqrmmBP'.
	Whit his binary you can create a RMQ structure. The parameters are documented in the code

References:
	Please, if you want to include this tool as part of your experiments, in your
	references include the paper [3] please.

[1]. J. Fischer and V. Heun. Space-efficient preprocessing schemes for range minimum 
queries on static arrays. SIAM '11.

[2]. K. Sadakane and G. Navarro. Fully-functional succinct trees. In Proceedings of 
the Twenty-First Annual ACM-SIAM Symposium on Discrete Algorithms, SODA '10.

[3]. H. Ferrada and G. Navarro. Improved Range Minimum Queries. In Proceedings of 
Data Compression Conference, DCC'16.
