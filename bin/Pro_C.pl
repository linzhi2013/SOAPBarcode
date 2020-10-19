#! usr/bin/perl -w
use strict;
use Getopt::Long;
=head1 Description
	
	--lib 	"lib type, f or s"
	--int 	"interval file of target region,formated in "site max min""
	--fas 	"sequence file in fasta format"
	--out 	"the filtered fasta"
	--help 	"print out this message"

=cut

my %mitocondon = (
        "TTT" => "F",   "TTC" => "F",   "TTA" => "L",
        "TTG" => "L",   "TCT" => "S",   "TCC" => "S",
        "TCA" => "S",   "TCG" => "S",   "TAT" => "Y",
        "TAC" => "Y",   "TAA" => "X",   "TAG" => "X",
        "TGT" => "C",   "TGC" => "C",   "TGA" => "W",
        "TGG" => "W",   "CTT" => "L",   "CTC" => "L",
        "CTA" => "L",   "CTG" => "L",   "CCT" => "P",
        "CCC" => "P",   "CCA" => "P",   "CCG" => "P",
        "CAT" => "H",   "CAC" => "H",   "CAA" => "Q",
        "CAG" => "Q",   "CGT" => "R",   "CGC" => "R",
        "CGA" => "R",   "CGG" => "R",   "ATT" => "I",
        "ATC" => "I",   "ATA" => "M",   "ATG" => "M",
        "ACT" => "T",   "ACC" => "T",   "ACA" => "T",
        "ACG" => "T",   "AAT" => "N",   "AAC" => "N",
        "AAA" => "K",   "AAG" => "K",   "AGT" => "S",
        "AGC" => "S",   "AGA" => "S",   "AGG" => "S",
        "GTT" => "V",   "GTC" => "V",   "GTA" => "V",
        "GTG" => "V",   "GCT" => "A",   "GCC" => "A",
        "GCA" => "A",   "GCG" => "A",   "GAT" => "D",
        "GAC" => "D",   "GAA" => "E",   "GAG" => "E",
        "GGT" => "G",   "GGC" => "G",   "GGA" => "G",
        "GGG" => "G",
);
my %stocondon = (
        "TTT" => "F",   "TTC" => "F",   "TTA" => "L",
        "TTG" => "L",   "TCT" => "S",   "TCC" => "S",
        "TCA" => "S",   "TCG" => "S",   "TAT" => "Y",
        "TAC" => "Y",   "TGT" => "C",   "TGC" => "C",   
        "TGA" => "W",   "GGG" => "G",   "GGA" => "G",
        "TGG" => "W",   "CTT" => "L",   "CTC" => "L",
        "CTA" => "L",   "CTG" => "L",   "CCT" => "P",
        "CCC" => "P",   "CCA" => "P",   "CCG" => "P",
        "CAT" => "H",   "CAC" => "H",   "CAA" => "Q",
        "CAG" => "Q",   "CGT" => "R",   "CGC" => "R",
        "CGA" => "R",   "CGG" => "R",   "ATT" => "I",
        "ATC" => "I",   "ATA" => "M",   "ATG" => "M",
        "ACT" => "T",   "ACC" => "T",   "ACA" => "T",
        "ACG" => "T",   "AAT" => "N",   "AAC" => "N",
        "AAA" => "K",   "AAG" => "K",   "AGT" => "S",
        "AGC" => "S",   "AGA" => "S",   "AGG" => "S",
        "GTT" => "V",   "GTC" => "V",   "GTA" => "V",
        "GTG" => "V",   "GCT" => "A",   "GCC" => "A",
        "GCA" => "A",   "GCG" => "A",   "GAT" => "D",
        "GAC" => "D",   "GAA" => "E",   "GAG" => "E",
        "GGT" => "G",   "GGC" => "G"
);
my %hydro=(
        "R"=>"-4.5",    "K"=>"-3.9",    "N"=>"-3.5",
        "D"=>"-3.5",    "Q"=>"-3.5",    "E"=>"-3.5",
        "H"=>"-3.2",    "P"=>"-1.6",    "Y"=>"-1.6",
        "W"=>"-0.9",    "S"=>"-0.8",    "T"=>"-0.7",
        "G"=>"-0.4",    "A"=>"1.8",     "M"=>"1.9",
        "C"=>"2.5",     "F"=>"2.8",     "L"=>"3.8",
        "V"=>"4.4",     "I"=>"4.5"
);
my ($Lib, $Fas, $Help, $Int, $Out);
GetOptions(
	"lib:s"=>\$Lib,
	"int:s"=>\$Int,
	"fas:s"=>\$Fas,
	"out:s"=>\$Out,
	"help:s"=>\$Help
);
die `pod2text $0` if ($Help || !defined $Lib || !defined $Fas || !defined $Out);
die "the hydrophily interval file is required for the full length lib" if (!defined $Int and ($Lib eq "f"));
my (%max,%min);
my $hal=0;
if (defined $Int){
	open INT, "$Int" or die $!;
	while(<INT>){
        	$hal++;
		chomp;
        	my @in=split /\t/;
        	$max{$in[0]}=$in[1];
        	$min{$in[0]}=$in[2];
	}
	close INT;
}
open FA, "$Fas" || die $!;
open OUT, ">$Out" || die $!;
if ($Lib eq "f"){
	$/="\>";
	<FA>;
	$/="\n";
	while(my $title=<FA>){
		chomp($title);
		$/="\>";
		chomp(my $seq=<FA>);
		$/="\n";
		$seq=~s/\n//g;
		my @seqm=split /NNN/,$seq;
		die "$title is not in the right form" unless (@seqm == 2);
		my $prof;
		my $lenf=length $seqm[0];
		for (my $i=1;$i<=$lenf-3;$i+=3){
			my $con=substr ($seqm[0],$i,3);
			if (exists $mitocondon{$con}){
				$prof.=$mitocondon{$con};
			}else{
				$prof.="X";
			}
		}
		my $pror;
		my $lenr=length $seqm[1];
		my $sta=($lenr % 3);
		for (my $i=$sta;$i<=$lenr-3;$i+=3){
			my $con=substr ($seqm[1],$i,3);
			if (exists $mitocondon{$con}){
				$pror.=$mitocondon{$con};
			}else{
				$pror.="X";
			}
		}
		next if (($prof=~/X/) or ($pror=~/X/));
		my $over=0;
		my $plenf=length $prof;
		my $plenr=length $pror;
		for (my $i=0; $i<=$plenf-3; $i++){
			my $subpro=substr($prof,$i,3);
			my $subval=value($subpro);
			if ($subval > $max{$i} || $subval < $min{$i}){
				$over++;
			}
		}
		my $rst=$hal-$plenr+2;
		for (my $i=0; $i<=$plenr-3;$i++){
			my $subpro=substr($pror,$i,3);
			my $subval=value($subpro);
			my $num=$rst+$i;
			if ($subval > $max{$num} || $subval < $min{$num}){
				$over++;
			}
		}
		if ($over <=1 ){
			print OUT ">$title\n$seq\n";
		}	
	}
}elsif($Lib eq "s"){
	$/="\>";
	<FA>;
	$/="\n";
	while(my $ti=<FA>){
		chomp($ti);
		$/="\>";
		chomp(my $seq=<FA>);
		$/="\n";
		$seq=~s/\n//g;
		my $seql=length $seq;
		my $seqtr=$seq;
		$seqtr=~tr/ATCGN/TAGCN/;
		$seqtr= reverse $seqtr;
		my $jud=0;
		TNT:for(my $start=0; $start<3; $start++) {
			my $couf=0;my $cour=0;my $juds=1;
			SEQ:for (my $i = $start; $i <= $seql - 3; $i += 3) {
				my $con = substr($seq, $i, 3);
				my $contr=substr($seqtr,$i,3);
#				if ($couf>=1 && $cour>=1){$juds=0;last SEQ;}
				if (!defined $stocondon{$con}){$couf++}
				if (!defined $stocondon{$contr}){$cour++}
				if ($couf>=1 && $cour>=1){$juds=0;last SEQ;}
			}
			if ($juds){
				$jud=1;
				last TNT;
			}
		}
		if ($jud){
			print OUT ">$ti\n$seq\n";
		}
	}


}else{
	die "the lib style is not in the right form, f for the full length lib and s for the shotgun lib\n";
}
sub value {
        my $strine=$_[0];
        my @b=split /\s*/,$strine;
        my $val=0;
        for my $num (0..$#b){
                if (exists $hydro{$b[$num]}){
                        $val+=$hydro{$b[$num]};
                }else{
                        print "$b[$num] is nothing make sense!\n";
                        $val+=0;
                }
        }
        return $val;
}
