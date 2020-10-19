#! usr/bin/perl -w
use strict;

use Getopt::Long;
=head1 Description

        --sam   "sam file"
        --fas   "sequence file in fasta format"
        --out   "the filtered fasta"
	--end 	"end length for chimera check | by default =100"
	--cut 	"cutoff value to remove abundance multiply variance | by default =1000"
        --hel 	"print out this message"

=cut
my ($Sam,$Fas,$Out,$Help,$End,$Cut);
GetOptions(
        "sam:s"=>\$Sam,
        "fas:s"=>\$Fas,
        "out:s"=>\$Out,
	"end:i"=>\$End,
	"cut:i"=>\$Cut,
        "help:s"=>\$Help
);
die `pod2text $0` if ($Help || !defined $Sam || !defined $Fas || !defined $Out);

$End ||= 100;
$Cut ||= 1000;
open IN, "$Sam" || die $!;
my %hash;
while (<IN>){
        chomp;
        if (/^\@SQ/){
                my @a=split;
                $a[1]=~s/SN://;
                $a[2]=~s/LN://;
                my @b;
                for (my $i=0; $i<=$a[2]-1; $i++){
                        $b[$i]=0;
                }
                $hash{$a[1]}= [@b];
        }
        next if (/^\@/);
        my @c=split;
        next if ($c[1] & 0x0004);
        my $len=length $c[9];
        for (my $ii=0; $ii<=$len-1; $ii++){
                my $as=$c[3]-1+$ii;
		if (defined $hash{$c[2]}[$as]){
                	$hash{$c[2]}[$as]++;
		}else{
			print "$c[2]\t$c[0]\n";
		}
        }
}
close IN;
open FAS, "$Fas" || die $!;
open OUT, ">$Out" || die $!;
$/="\>";
while(<FAS>){
	chomp;
	next if ($_ eq "");
	my $title=(split /\s+/)[0];
	my $seq=(split /\n/,$_,2)[1];
	$seq=~s/\n//g;
	next unless (exists $hash{$title});
	my $tot=0; my $fab=0; my $rab=0;
	for my $i (0..$#{$hash{$title}}){
		$tot+=$hash{$title}[$i];
		$fab+=$hash{$title}[$i] if ($i<$End);
		$rab+=$hash{$title}[$i] if ($i>($#{$hash{$title}}-$End));
	}
	next if ($tot==0 or $fab==0 or $rab==0);
	my $lenth=$#{$hash{$title}}+1;
	my $abun=int($tot/$lenth);
	my $fav=($fab/$End);
	my $rav=($rab/$End);
	next if (($fav/$rav >= $Cut) || ($rav/$fav >= $Cut));
	print OUT ">$title;size=$abun\n$seq\n"
}

close FAS;
close OUT;

print "all done!\n";
