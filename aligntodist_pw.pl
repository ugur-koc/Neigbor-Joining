#!/usr/bin/perl

# This script creates distance matrixes using pairwise global alignment scores
# one may need to modify diroctries
# input: none
# output: it will print out the ditance matrix that is ready to pass our NJ implementation
# @auhtor Ugur Koc
# @date Dec.8.2015

use strict;
use warnings;
use Math::BigFloat;

my $workDir = "/Users/ugurmeryem/Dropbox/701CMSC/Neigbor-Joining/data";

my @seqnames;
open F_NAMES,"<$workDir/85VASTdomains_names.txt" or die $!;
my $line = <F_NAMES>;
@seqnames = split(",", $line);
close F_NAMES;

my @alignmentFiles = ();
opendir(DIR, "$workDir/g_alignments") or die $!;
while (my $file = readdir(DIR)) {
      push(@alignmentFiles, $file);
}

my %score;
my $seq1;
my $seq2;
foreach my $file (@alignmentFiles){
   open ALIGN,"<$workDir/g_alignments/$file" or die $!;
   while (my $line = <ALIGN>){
      if($line =~ / 1: /){
         my @tmp = split(": ", $line);
         chomp @tmp;
         $seq1 = $tmp[1];
         #print "seq1: $seq1";
      }elsif($line =~ / 2: /){
         my @tmp = split(": ", $line);
         chomp @tmp;
         $seq2 = $tmp[1];
         #print "seq2: $seq2";
      }elsif($line =~ / Score: /){
         my @tmp = split(": ", $line);
         chomp @tmp;
         $score{$seq1}{$seq2}=$tmp[1];
         #print "score: $score{$seq1}{$seq2}";
         last;
      }
   }
}
#exit;
print "   85\n";
print "@seqnames\t";
print "\n";
foreach my $seq1 (@seqnames) {
   foreach my $seq2 (@seqnames) {
      if($seq1 eq $seq2) {
         print "0\t";
      } else {
         my $distance = 1000/$score{$seq1}{$seq2};
         print "$distance\t";
      }
   }
   print "\n";
}