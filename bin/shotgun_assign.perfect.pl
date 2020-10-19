#! usr/bin/perl -w
use strict;
use Getopt::Long;

=head1 Description

	Usage: perl 

	--pri 	"primer set with index ahead"
	--fas 	"input file as fasta format"
	--out 	"the prefix of all the output file"
	--mis 	"the maximun mismatch allowed for the primer | by default =0"
	--len	"the length of the index | by default =0, set as 0 when no index"
	--help 	"print out this message"

=cut

my%primer=(
        "V" => "(A|C|G)",
        "D" => "(A|T|G)",
        "B" => "(T|G|C)",
        "H" => "(A|T|C)",
        "W" => "(A|T)",
        "S" => "(C|G)",
        "K" => "(T|G)",
        "M" => "(A|C)",
        "Y" => "(C|T)",
        "R" => "(A|G)",
        "N" => "(A|T|C|G)",
);
my ($Pri, $Fas, $Out, $Mis, $Len, $Help);

GetOptions(
	"pri:s"=>\$Pri,
	"fas:s"=>\$Fas,
	"out:s"=>\$Out,
	"mis:i"=>\$Mis,
	"len:i"=>\$Len,
	"help"=>\$Help
);

die `pod2text $0` if ($Help || !defined $Pri || !defined $Fas || !defined $Out);
$Mis ||= 0;
$Len ||= 0;
die "mis value should be 0 which leave no space for argument" unless ($Mis == 0);
my $priwc;
$priwc = `less -S $Pri|wc -l`;
die "the primer set need to have each forward and reverse primer in demanded format"
unless ($priwc%4 == 0);

die "the primer file contain more than 2 primer sets, however, without index length hava been seted" if ($priwc > 4 && $Len==0);
my (%prif,%prir);
open PRI, "$Pri" || die $!;

my ($plenf, $plenr);

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
                my $forp= convert ($for);
                $prif{$ind1}=$forp;
#		my @a=split /\s*/,$for;
#		for (my $i=0; $i<=$#a; $i++){
#			if (exists($primer{$a[$i]})){
#				push @{$prif{$ind1}}, $primer{$a[$i]}
#			}else{
#				push @{$prif{$ind1}}, $a[$i]
#			}
#		}
                my $revp= convert ($rev);
                $prir{$ind2}=$revp;
#		my @b=split /\s*/, $rev;
#		for (my $i=0; $i<=$#b; $i++){
#			if (exists($primer{$b[$i]})){
#				push @{$prir{$ind2}}, $primer{$b[$i]}
#			}else{
#				push @{$prir{$ind2}}, $b[$i]
#			}
#		}
	}else{
		my $ind=1;
                my $forp= convert ($for);
                $prif{$ind}=$forp;
                my $revp= convert ($rev);
                $prir{$ind}=$revp;
#		my @a=split /\s*/,$for;
#               for (my $i=0; $i<=$#a; $i++){
#                        if (exists($primer{$a[$i]})){
#                                push @{$prif{$ind}}, $primer{$a[$i]}
#                        }else{
#                                push @{$prif{$ind}}, $a[$i]
#                        }
#                }
#                my @b=split /\s*/, $rev;
#                for (my $i=0; $i<=$#b; $i++){
#                        if (exists($primer{$b[$i]})){
#                                push @{$prir{$ind}}, $primer{$b[$i]}
#                        }else{
#                                push @{$prir{$ind}}, $b[$i]
#                        }
#                }
	}
}
sub convert {
        my $old = $_[0];
        my @a = split /\s*/,$old;
        my $primer;
        for (my $i=0; $i<=$#a; $i++){
                if (exists $primer{$a[$i]}){
                        $primer.=$primer{$a[$i]};
                }else{
                        $primer.=$a[$i];
                }
        }
        return $primer;
}
#my $subc;
#if ($plenr>=$plenf){
#	$subc=$plenr
#}else{
#	$subc=$plenf
#}
#my $fcut=$plenf-$Mis;
#my $rcut=$plenr-$Mis;
open FA, "$Fas" || die $!;
my %FH;
open MIL, ">$Out\_list" || die $!;
for my $key (keys %prif){
	open ($FH{"$key"}, ">", "$Out\_$key\.fasta") || die $!;
	print MIL "$Out\_$key\.fasta\n";
	}
close MIL;
TNT:while (<FA>){
	chomp;
	my $ti=(split /\s+/)[0];
	chomp (my $seq=<FA>);
	my $seqt=$seq;
	$seqt=~tr/ATCGN/TAGCN/;
	$seqt= reverse $seqt;
#	my $seq_s=substr($seq,0,$subc);
#	my $seqt_s=substr($seqt,0,$subc);
#	my @ss=split /\s*/,$seq_s;
#	my @ts=split /\s*/,$seqt_s;
	my $noj;
	QEC:for my $key (keys %prif){
		my $ffj=0;my $frj=0;my $rfj=0;my $rrj=0;
		$noj=0;
                if ($seq=~s/^$prif{$key}//){
                	$ffj=1;
                }
		if($seq=~s/^$prir{$key}//){
			$frj=1;
		}
                if ($seqt=~s/^$prif{$key}//){
			$rfj=1;
		}
		if ($seqt=~s/^$prir{$key}//){
			$rrj=1;
		}
#		if (($ffj==1 && $frj>=$fcut) or ($ffj>=$fcut && $rfj >= $rcut) or ($frj>=$fcut && $rfj >= $rcut) or ($rfj >= $rcut && $rrj >= $rcut) or ($ffj>=$fcut && $rrj >= $rcut) or ($frj>=$fcut&& $rrj >= $rcut)){next TNT};
		my $delj=$ffj+$frj+$rfj+$rrj;
		next TNT if ($delj >= 2);
		if ($ffj==1){
			print {$FH{$key}} "$ti\_f\n$seq\n";
			last QEC;
		}elsif($frj==1){
			print {$FH{$key}} "$ti\_r\n$seq\n";
			last QEC;
		}elsif ($rfj ==1){
			print {$FH{$key}} "$ti\_f\n$seqt\n";
			last QEC;
		}elsif($rrj ==1){
			print {$FH{$key}} "$ti\_r\n$seqt\n";
			last QEC;
		}else{
			$noj=1;
		}
	}
	if ($noj){
		for my $key (keys %prif){
			print {$FH{$key}} "$ti\n$seq\n";
		}
	}
}

