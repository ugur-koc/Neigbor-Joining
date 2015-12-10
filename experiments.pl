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

my $data="85VAST";
my $groundTruth="data/$data-groundtruth.tree";

my $exp=$ARGV[0];
if (not defined $exp) {
   $exp=1;
}
if($exp==1){
   my @protMethods = ("Uncorrected", "Jukes-Cantor", "KimuraProtein");
   my @gapScores = ("1", "25", "75");
   my $size=scalar @gapScores;
   for(my $i=0; $i < 3; $i++){
      if($i==3){ $size=1;}
      for(my $j=0; $j < $size; $j++){
         
         print "method:$protMethods[$i],gapScore:$gapScores[$j],";
         #Following like produces distance matrix ready to feed into our NJ implementation
         my $distanceFile="data/$data-$i-$j.dist";
         my $start = sprintf("%.3f", Time::HiRes::time());
         `./aligntodist.pl $data.aln tmp $i $gapScores[$j] > $distanceFile`;
         my $time = sprintf("%.3f", Time::HiRes::time() - $start);
         print "Tdist:$time,";
         
         #Following produces the tree given the distance matrix; tree is in $nj_out file
         my $nj_out="data/$data-tree-$i-$j.tmp";
         $start = sprintf("%.3f", Time::HiRes::time());
         `./main $distanceFile $nj_out`;
         $time = sprintf("%.3f", Time::HiRes::time() - $start);
         print "Ttree:$time,";
         
         #Following line should help me compare trees; visualize them produce some score accuracy etc.
         my $ct_out="data/$data-$i-$j.tree";
         `./compareTrees.py $groundTruth $nj_out > $ct_out`;
         
         open F_TREE,"<$ct_out" or die $!;
         my $line = <F_TREE>;
         print "$line";
         close F_TREE;
      }
   }
}elsif($exp==2){
   my $distanceFile="data/$data-pw.dist";
   my $start = sprintf("%.5f", Time::HiRes::time());
   `./aligntodist_pw.pl > $distanceFile 2>&1`;
   my $time = sprintf("%.3f", Time::HiRes::time() - $start);
   print "method:pairwise,gapScore:10,Tdist:$time,";
   
   my $nj_out="data/$data-tree-pw.tmp";
   $start = sprintf("%.3f", Time::HiRes::time());
   `./main $distanceFile $nj_out`;
   $time = sprintf("%.3f", Time::HiRes::time() - $start);
   print "Ttree:$time,";
   
   my $ct_out="data/$data-pw.tree";
   `./compareTrees.py $groundTruth $nj_out > $ct_out`;
   
   open F_TREE,"<$ct_out" or die $!;
   my $line = <F_TREE>;
   print "$line";
   close F_TREE;
}elsif($exp==3){
   
}
`rm -rf data/*.tmp`;