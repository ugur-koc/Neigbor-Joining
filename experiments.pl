#!/usr/bin/perl

use strict;
use warnings;
use Time::HiRes;
use Math::BigFloat;

my @protMethods = ("0", "1", "2");
my @gapScores = ("0", "0.1", "1");

my $data="85VASTdomains_clean.aln";
for(my $i=0; $i < 3; $i++){
   for(my $j=0; $j < 3; $j++){
      my $start = Time::HiRes::time();
      `./aligntodist.pl $data tmp $protMethods[$i] $gapScores[$j] > $data-$i-$j.txt`;
      my $time = Time::HiRes::time() - $start;
      print "$i, $j : $time\n"
   }
}