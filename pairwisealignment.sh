#/bin/sh

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