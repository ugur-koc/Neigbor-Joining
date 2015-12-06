#!/usr/bin/perl

use strict;
use warnings;
use Time::HiRes;
use Math::BigFloat;

my $workDir = "/Users/ugurmeryem/Dropbox/701CMSC/Neigbor-Joining/data/"; #$ARGV[0];
my $data = "$workDir/85VASTdomains.aln";
my $distmatOut = "$workDir/tmp.txt";
my $outfile = '$workDir/$data.dm';

`distmat -sequence $data -protmethod 1 -outfile $distmatOut`;

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