CC = gcc
CP = python3
LIBS = -lgsl -lgslcblas -lm

master:Main.c
	$(CC) -o main.o Main.c $(LIBS)
	./main.o
	#$(CP) Plot.py

	@echo ""
	@echo "<-------------------------	.THE END.	------------------------->"
	@echo ""

trajdata:TrajFiles.py
	$(CP) TrajFiles.py

test:binterpole_tst.c
	$(CC) -o binterpole_tst.o binterpole_tst.c $(LIBS)
	./binterpole_tst.o

bin:RetriveBinData.c
	$(CC) -o bin.o RetriveBinData.c
	./bin.o

main:Main.c 
	$(CC) -o main.o Main.c $(LIBS)
	./main.o

abs:get_L_abs_vol.c
	$(CC) -o abs.o get_L_abs_vol.c $(LIBS)
	./abs.o

