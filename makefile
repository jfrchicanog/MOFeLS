CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio1262/cplex
INCLUDE=$(CPLEX_HOME)/include
STATIC_LIB=$(CPLEX_HOME)/lib/x86-64_linux/static_pic
DYNAMIC_LIB=$(CPLEX_HOME)/bin/x86-64_linux
CPLEX_LIB_STATIC=cplex
CPLEX_LIB_DYNAMIC=cplex1262
LFLAGS=-lrt -lm -lpthread
CFLAGS=-DIL_STD -I$(INCLUDE) -O0 -g3 -c -fmessage-length=0 -std=c++0x
CC=g++

MOFeLS:	main_MOFeLS.o
	$(CC) -L$(STATIC_LIB) -o MOFeLS ./main_MOFeLS.o -l$(CPLEX_LIB_STATIC) $(LFLAGS)

main_MOFeLS.o:	main_MOFeLS.cpp
	$(CC) $(CFLAGS) main_MOFeLS.cpp

clean:
	rm main_MOFeLS.o MOFeLS
