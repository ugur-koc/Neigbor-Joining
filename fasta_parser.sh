#!/bin/sh

# This script parses .fasta file that contains multible entries
# and creates a .fata file for each sequence
# input: input .fasta file (one may need to modify diroctries in lines 15 and 16)
# output: .fasta file for each sequence
# @auhtor Ugur Koc
# @date Dec.8.2015

counter=1;
while IFS='' read -r line || [[ -n "$line" ]]; do
	rrIN=(${line//	/ })
    seqid=${rrIN[0]};
    seq=${rrIN[1]};
    echo "$seqid" > data/fastas/85VASTdomains-$counter.fasta
    echo "$seq" >> data/fastas/85VASTdomains-$counter.fasta
    counter=$((counter+1));
done < "$1"