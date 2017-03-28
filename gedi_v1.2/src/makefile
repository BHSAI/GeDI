#
CC=g++
GSLDIR=/home/hwoo/local
CFLAGS=-I$(GSLDIR)/include -c -Wall -g
LDFLAGS=-L$(GSLDIR)/lib -lgsl -lgslcblas -g
#CFLAGS=-c -Wall -g
#LDFLAGS=-lgsl -lgslcblas -g

all: gedi

gedi: main.o il.o il_lr.o cl.o cl_lr.o qt.o qt_lr.o
	$(CC) main.o il.o il_lr.o cl.o cl_lr.o qt.o qt_lr.o $(LDFLAGS) -o gedi

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

qt.o: qt.cpp
	$(CC) $(CFLAGS) qt.cpp

qt_lr.o: qt_lr.cpp
	$(CC) $(CFLAGS) qt_lr.cpp

clean: 
	rm -rf *.o gedi
