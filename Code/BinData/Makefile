CC = gcc
CP = python3
LIBS = -lgsl -lgslcblas -lm

bin:RetriveBinData.c
	$(CC) -o bin.o RetriveBinData.c
	./bin.o
	rm -r *.o

