#! usr/bin/perl -w
use strict;
=head1 Description

        Usage: perl extract.pl <parameter>

        --fq1    	"xx.1.fq"
        --fq2    	"xx.2.fq"
        --pri		"primer set with index ahead"
        --out		"prefix of the output file"
        --len     	"the length of index | by default =0"
	--mis		"the allowed mismatch of index+primer | by default =1"
	--qua		"the maximun B qualitiy base allowed | by default =5"
	--help   	"print out this information"

=cut

use Getopt::Long;
use FindBin qw($Bin $Script);
use strict;
use PerlIO::gzip;

my ($Fq1,$Fq2,$Pri,$Out,$Len,$Mis,$Qua,$Help);

my %primer=(
        "V" => "ACG",
        "D" => "ATG",
        "B" => "TGC",
        "H" => "ATC",
        "W" => "AT",
        "S" => "CG",
        "K" => "TG",
        "M" => "AC",
        "Y" => "CT",
        "R" => "AG",
        "N" => "ATCG",
       );


GetOptions(
        "fq1:s"=>\$Fq1,
        "fq2:s"=>\$Fq2,
        "out:s"=>\$Out,
	"pri:s"=>\$Pri,
        "len:i"=>\$Len,
        "mis:i"=>\$Mis,
	"qua:i"=>\$Qua,
        "help"=>\$Help
);

die `pod2text $0` if ($Help || !defined ($Fq1) || !defined ($Fq2)||!defined ($Pri));
$Len ||= 0;
$Mis=1 if (!defined $Mis);
$Qua=5 if (!defined $Qua);

if ($Fq1=~/gz$/){
                open (FFQ, "<:gzip",$Fq1) || die $!;
        }else{
                open (FFQ, $Fq1) || die $!;
        }
if ($Fq2=~/gz$/){
                open (RFQ,"<:gzip",$Fq2) || die $!;
        }else{
                open (RFQ, $Fq2) || die $!;
        }

my $priwc;
$priwc = `less -S $Pri|wc -l`;
die "the primer file need to have each forward and reverse primer"
unless ($priwc%4 == 0);

die "the primer file contain more than 2 primer sets, however, without index length hava been seted" if ($priwc > 4 && $Len==0);
print "the len is $Len\n";
open PRI, "$Pri" || die $!;
my ($plenf, $plenr,%prif,%prir);
while(<PRI>){
        chomp;
        chomp(my $for=<PRI>);
        $plenf=length $for;
        <PRI>;
        chomp(my $rev=<PRI>);
        $plenr=length $rev;
        if ($Len != 0){
                my $ind1=substr($for,0,$Len);
                my $ind2=substr($rev,0,$Len);
                die "the index of the forwad and reverse primers are not identical or have index duplication between differnet primers" if (($ind1 ne $ind2) or exists $prif{$ind1});
                my @a=split /\s*/,$for;
                for (my $i=0; $i<=$#a; $i++){
                        if (exists($primer{$a[$i]})){
                                push @{$prif{$ind1}}, $primer{$a[$i]}
                        }else{
                                push @{$prif{$ind1}}, $a[$i]
                        }
                }
                my @b=split /\s*/, $rev;
                for (my $i=0; $i<=$#b; $i++){
                        if (exists($primer{$b[$i]})){
                                push @{$prir{$ind2}}, $primer{$b[$i]}
                        }else{
                                push @{$prir{$ind2}}, $b[$i]
                        }
                }
        }else{
                my $ind=1;
                my @a=split /\s*/,$for;
                for (my $i=0; $i<=$#a; $i++){
                        if (exists($primer{$a[$i]})){
                                push @{$prif{$ind}}, $primer{$a[$i]}
                        }else{
                                push @{$prif{$ind}}, $a[$i]
                        }
                }
                my @b=split /\s*/, $rev;
                for (my $i=0; $i<=$#b; $i++){
                        if (exists($primer{$b[$i]})){
                                push @{$prir{$ind}}, $primer{$b[$i]}
                        }else{
                                push @{$prir{$ind}}, $b[$i]
                        }
                }
        }
}
close PRI;
my $fcut=$plenf-$Mis;
my $rcut=$plenr-$Mis;
my %FH;
open FUL, ">$Out\_list" || die $!;
for my $key (keys %prif){
        open ($FH{"$key"}, ">", "$Out\_$key\.fasta") || die $!;
        print FUL "$Out\_$key\.fasta\n";
}
close FUL;
my $subc;
if ($plenr>=$plenf){
        $subc=$plenr
}else{
        $subc=$plenf
}
my $seqnum=0;
TNT:while(<FFQ>){
        $seqnum++;
	chomp;
        chomp(my $seqf=<FFQ>);
        <FFQ>;
        chomp(my $quaf=<FFQ>);
        <RFQ>;
        chomp (my $seqr=<RFQ>);
        <RFQ>;
        chomp(my $quar=<RFQ>);
	next TNT if ($seqf=~/NN/ or $seqr=~/NN/);
	my $temlenf=length $seqf;
	my $temlenr=length $seqr;
	next TNT unless ($temlenf >= 100 && $temlenr >= 100);
	my $Bnumf=$quaf=~s/B//g;
        my $Bnumr=$quar=~s/B//g;
        if ($Bnumf<= $Qua && $Bnumr<=$Qua){
		my $seqf_s=substr($seqf,0,$subc);
		my $seqr_s=substr($seqr,0,$subc);
		my @sfs=split /\s*/, $seqf_s;
		my @srs=split /\s*/, $seqr_s;
		QEC:for my $key (keys %prif){
			my $ffj=0;my $frj=0;my $rfj=0;my $rrj=0;
			for (my $i=0; $i<= $#{$prif{$key}}; $i++){
				if($prif{$key}[$i]=~m/$sfs[$i]/){
                        	        $ffj++;
	                        }
        	                if($prif{$key}[$i]=~m/$srs[$i]/){
                	                $frj++;
                        	}
			}
			for (my $i=0; $i<= $#{$prir{$key}}; $i++){
                        	if ($prir{$key}[$i]=~m/$sfs[$i]/){
                                	$rfj++;
                        	}
                        	if ($prir{$key}[$i]=~m/$srs[$i]/){
                                	$rrj++;
                        	}
                	}
			my $seqm;
			next TNT if (($ffj>=$fcut && $frj>=$fcut) or ($rfj>=$rcut && $rrj>=$rcut));
			if ($ffj>=$fcut && $rrj>=$rcut){
				my $seqfw=substr($seqf,$plenf);
				my $seqrw=substr($seqr,$plenr);
				$seqrw = reverse $seqrw;
				$seqrw=~tr/ATCGN/TAGCN/;
				$seqm = "$seqfw"."NNN"."$seqrw";
#				print "f\t$ffj\t$fcut\t$rrj\t$rcut\t$_\n";
				print {$FH{$key}} ">$key\_$seqnum\n$seqm\n";				
			}elsif($rfj>=$rcut && $frj>=$fcut){
				my $seqfw=substr($seqf,$plenr);
				my $seqrw=substr($seqr,$plenf);
				$seqfw = reverse $seqfw;
				$seqfw=~tr/ATCGN/TAGCN/;
				$seqm = "$seqrw"."NNN"."$seqfw";
#				print "r\t$rfj\t$rcut\t$frj\t$fcut\t$_\n";
				print {$FH{$key}} ">$key\_$seqnum\n$seqm\n";
			}
		}
	}
}
close FFQ;
close RFQ;
#print "assign has been done\n";

#open LIS, "$Out\_flist" || die $!;
#my $binpath=$Bin;
#while(my $file=<LIS>){
#	chomp ($file);
#	my $outd="$file"."dup";
#	`$binpath/bin/usearch7 -derep_fulllength $file -output $outd -sizeout -strand both -minuniquesize 2`;
#}
#close LIS;
print "all done\n";
