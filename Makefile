CC=g++
CFLAGS=-O9 -DNDEBUG

all: index

Basic_drf64.o: includes/Basic_drf64.cpp
	$(CC) $(CFLAGS) -c includes/Basic_drf64.cpp

Chameleon.o: Chameleon.cpp
	$(CC) $(CFLAGS) -c Chameleon.cpp

ConfigFile.o: ConfigFile.cpp
	$(CC) $(CFLAGS) -c ConfigFile.cpp

RangeMMTree64.o: RangeMMTree64.cpp
	$(CC) $(CFLAGS) -c RangeMMTree64.cpp

RMQRMM64.o: RMQRMM64.cpp
	$(CC) $(CFLAGS) -c RMQRMM64.cpp

index: Basic_drf64.o Chameleon.o ConfigFile.o RangeMMTree64.o RMQRMM64.o
	ar rc drflib64.a Basic_drf64.o Chameleon.o ConfigFile.o  RangeMMTree64.o RMQRMM64.o

test: 
	@$(CC) $(CFLAGS) test.cpp -o myTest drflib64.a 

clean:
	-rm *~ *.o *.bak 
cleanall:
	-rm *~ *.o *.bak .a
