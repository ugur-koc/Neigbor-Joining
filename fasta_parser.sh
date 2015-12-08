#!/bin/sh

counter=1;
while IFS='' read -r line || [[ -n "$line" ]]; do
	rrIN=(${line//	/ })
    seqid=${rrIN[0]};
    seq=${rrIN[1]};
    echo "$seqid" > data/fastas/85VASTdomains-$counter.fasta
    echo "$seq" >> data/fastas/85VASTdomains-$counter.fasta
    counter=$((counter+1));
done < "$1"