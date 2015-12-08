#!/usr/bin/perl

# This script runs the experiment conducted in this project
# Run this script in the same directory with main file
# input: none (one should modify the directory variables accordingly)
# output: it will print out following for each construction approach;
#           time to compute distance matrix,
#           time to compute phylogenetic tree,
#           and an accuracy score
# @auhtor Ugur Koc
# @date Dec.8.2015

use strict;
use warnings;
use Time::HiRes;
use Math::BigFloat;

my @protMethods = ("0", "1", "2");
my @gapScores = ("0", "0.1", "1");

my $data="85VASTdomains_clean";
my $groundTruth="85VASTdomains.newick";
for(my $i=0; $i < 3; $i++){
   for(my $j=0; $j < 3; $j++){
      my $start = Time::HiRes::time();
      my $distanceFile="data/$data-dist-$i-$j.txt";
      
      #Following like produces distance matrix ready to feed into our NJ implementation
      `./aligntodist.pl $data.aln tmp $protMethods[$i] $gapScores[$j] > $distanceFile`;
      my $time = Time::HiRes::time() - $start;
      print "DM computation time for method:$i, and gap score:$j --> $time\n"
      my $nj_out="$data-tree-$i-$j.txt";
      
      #Following produces the tree given the distance matrix; tree is in $nj_out file
      `./main $distanceFile $nj_out`;
      
      #Following line should help me compare trees; visualize them produce some score accuracy etc.
      `python2.7 compareTrees.py $groundTruth $nj_out`;
   }
}