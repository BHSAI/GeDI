#
CC=g++
#CC=icpc
#CFLAGS=-c -Wall -g
#CFLAGS=-c -Wall -DMPI  # for MPI
GSLDIR=/home/hwoo/local
CFLAGS=-I$(GSLDIR)/include -c -Wall -g
#CFLAGS=-I$(GSLDIR)/include -c -Wall -O2
#LDFLAGS=-L$(GSLDIR)/lib -lgsl -lgslcblas -O2 -lmpi -lmpi++
LDFLAGS=-L$(GSLDIR)/lib -lgsl -lgslcblas -g
#LDFLAGS=-L$(GSLDIR)/lib -lgsl -lgslcblas -O2

all: gedi

gedi: main.o il.o il_lr.o cl.o cl_lr.o
	$(CC) main.o il.o il_lr.o cl.o cl_lr.o $(LDFLAGS) -o gedi

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

il.o: il.cpp
	$(CC) $(CFLAGS) il.cpp

il_lr.o: il_lr.cpp
	$(CC) $(CFLAGS) il_lr.cpp

cl.o: cl.cpp
	$(CC) $(CFLAGS) cl.cpp

cl_lr.o: cl_lr.cpp
	$(CC) $(CFLAGS) cl_lr.cpp

clean: 
	rm -rf *.o gedi
