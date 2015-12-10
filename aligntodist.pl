#!/usr/bin/perl

# This script creates distance matrixes using distmat program from emboss package
# one may need to modify diroctries
# input: .fasta file that contains the multible alignment
#        output file name to pass to distmat
#        0, 1, or 2 to select the model of evolution
#        gap score
# output: it will print out the ditance matrix that is ready to pass our NJ implementation
# @auhtor Ugur Koc
# @date Dec.8.2015

use strict;
use warnings;
use Time::HiRes;
use Math::BigFloat;

my $workDir = "/Users/ugurmeryem/Dropbox/701CMSC/Neigbor-Joining/data";
my $data = "$workDir/$ARGV[0]";
my $distmatOut = "$workDir/$ARGV[1].tmp";
my $protMethod = $ARGV[2];
my $gapScore = $ARGV[3];

`touch $distmatOut`;
`distmat -sequence $data -protmethod $protMethod -gapweight $gapScore -outfile $distmatOut > data/z.tmp 2>&1`;

my %myMap = ();
open FILE,"<$distmatOut" or die $!;
my @distanceMatrix = ();
my $counter = 0;
my %seqIDMap = ();
while (my $line = <FILE>){
   $counter = $counter + 1;
   if($counter > 7){
      my $lasttab = rindex($line, '	');
      my @scores = split(' ', substr($line, 0, $lasttab));
      chomp @scores;
      my @seqName = split(' ', substr($line, $lasttab));
      $seqIDMap{$seqName[1]}=$seqName[0];
      $myMap{$seqName[0]} = \@scores;
   }
}

print "   85\n";
my $size = keys %seqIDMap;
for(my $i=1; $i <= $size; $i++){
   print $seqIDMap{$i}."\t";
}
print "\n";
my $content = "";
for(my $i=1; $i <= $size; $i++){
   my $scoresRef = $myMap{$seqIDMap{$i}};
   my $len = scalar @{ $scoresRef };
   for(my $j=1; $j <= $size; $j++){
      if ($j-1 < $size-$len){
         my $ugur = $myMap{$seqIDMap{$j}};
         $content = $content."".$ugur->[$i-$j]."\t";
      }else{
         $content = $content."".$scoresRef->[$j-1-$size+$len]."\t";
      }
   }
   $content = $content.""."\n";
}
$content=~s/nan/99999/g;
print $content;