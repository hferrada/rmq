CC=g++
CFLAGS=-O9 -DNDEBUG

all: index

Basic_rmq.o: includes/Basic_rmq.cpp
	$(CC) $(CFLAGS) -c includes/Basic_rmq.cpp

RMQRMM64.o: RMQRMM64.cpp
	$(CC) $(CFLAGS) -c RMQRMM64.cpp

index: Basic_rmq.o RMQRMM64.o
	ar rc rmqrmm.a Basic_rmq.o RMQRMM64.o

test: 
	@$(CC) $(CFLAGS) test.cpp -o myRMQ rmqrmm.a 

clean:
	-rm *~ *.o *.bak 
cleanall:
	-rm *~ *.o *.bak .a