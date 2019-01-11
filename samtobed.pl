#!/usr/bin/perl
###############################################################
# Program Name: sam_to_bed.pl             Version: 1.0
#
# Description: easily convert SAM file format to BED format
#
#
# Author: Shan Sabri          Date: 2017-11-21
#
#
# Revision History
#
# Version                 Date                  Who
#-------------------------------------------------------------
#  INIT                 2017-11-21             Shan Sabri  
###############################################################
#
# 
use strict;
use warnings;

my $f = $ARGV[0];

open IN, $f || die "Can't open the file!";
my $k=0;
while(<IN>){
	chomp;
	my @line = split /\t/;
	my $str = "-";
	$str = "+" if ($line[1] == 0); 
	my @ll = split (/\,/,$line[11]);
	foreach my $l (@ll){
		my ($s,$e) = split (/\:/,$l);
		--$s;
		print "$line[2]\t$s\t$e\t$line[0]\t0\t$str\n";
		#if ($k%2 == 0){
		#print "$line[2]\t$s\t$e\t$line[0]^1\t0\t$str\n";
		#}else{
		#print "$line[2]\t$s\t$e\t$line[0]^2\t0\t$str\n";
		#}
	}
	++$k;
}
close IN;
