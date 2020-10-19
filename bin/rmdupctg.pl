#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV!=2) {
	print "Usage:perl $0 [in] [out]\n";
	exit;
}

open (I,"$ARGV[0]") || die "can't open file $ARGV[0]\n";
open (O,">$ARGV[1]") || die "can't open file $ARGV[1]\n";
$/=">";
my %hash;
while (<I>) {
	chomp;
	next if ($_ eq "");
	my $id=(split /\s+/,$_)[0];
	my $seq=(split (/\n/,$_,2))[1];
	$seq=~s/\n//g;
	my $seqt=$seq;
	$seqt=~tr/ATCGN/TAGCN/;
	$seqt=reverse $seqt;
	if ((exists $hash{$seq}) or (exists $hash{$seqt})) {
		next;
	}else{
		$hash{$seq}=$id;
	}
}

foreach my $i (keys %hash) {
	print O ">$hash{$i}\n$i\n";
}
