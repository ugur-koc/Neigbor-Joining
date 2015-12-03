all:
	g++ main.cpp joiner.cpp -std=c++0x -O2 -Wall -o main


distance:
	gcc -c -g -fprofile-arcs -ftest-coverage phylip.c
	gcc -c -g -fprofile-arcs -ftest-coverage protdist.c
	gcc -g -fprofile-arcs -ftest-coverage protdist.o phylip.o -lm -o protdist

clean:
	rm -f *.o *.gcno *.gcda