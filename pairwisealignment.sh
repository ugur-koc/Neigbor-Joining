#/bin/sh

# This script creates pair alignments
# one may need to modify directries and run it in the same directory with main file
# input: none (just check the directorys used in the forr loops)
# output: it will creat an alinment file for each pair of sequences
# @auhtor Ugur Koc
# @date Dec.8.2015

counter1=0;
for i in data/fastas/85VASTdomains-*.fasta; 
do 
	counter2=0;
	counter1=$((counter1+1));
	for j in data/fastas/85VASTdomains-*.fasta; 
	do 
	    counter2=$((counter2+1));
		needle $i $j stdout -gapopen 10.0 -gapextend 0.5 >> data/distances/85VASTdomains-$counter1-$counter2.dist; 
	done; 
done;