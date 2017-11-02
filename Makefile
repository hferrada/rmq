CC=g++
CFLAGS=-std=c++11 -O3 -DNDEBUG -march=native

all: index

Basic_rmq.o: includes/Basic_rmq.cpp
	$(CC) $(CFLAGS) -c includes/Basic_rmq.cpp

RMQRMM64.o: RMQRMM64.cpp
	$(CC) $(CFLAGS) -c RMQRMM64.cpp

index: Basic_rmq.o RMQRMM64.o
	ar rc rmqrmmBP.a Basic_rmq.o RMQRMM64.o

test: 
	@$(CC) $(CFLAGS) rmqrmmBP.cpp -o rmqrmmBP rmqrmmBP.a 

clean:
	-rm *~ *.o *.bak 
cleanall:
	-rm *~ *.o *.bak .a
