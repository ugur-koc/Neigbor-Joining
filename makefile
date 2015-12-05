all:
	g++ main.cpp joiner.cpp -std=c++0x -O2 -Wall -o main

distance:
	gcc -c phylip.c
	gcc -c protdist.c
	gcc protdist.o phylip.o -lm -o protdist

clean:
	rm -f *.o protdist main
