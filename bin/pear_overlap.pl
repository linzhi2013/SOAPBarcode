#! usr/bin/perl -w
use strict;
use FindBin qw($Bin $Script);
use Getopt::Long;

=head1 description

	perl $0 <parameters>
	
	-for 	forward 1.fq
	-rev 	reverse 2.fq
	-out 	the output
	-min 	the minimum length of the assembled PE
	-mio 	the minimum of overlaped legnth
	-cpt 	the number of threads
	-help 	print out this message

=cut

my ($For, $Rev, $Mao, $Min, $Mio, $Out, $Cpt,$Help);

GetOptions(
	"for:s"=>\$For,
	"rev:s"=>\$Rev,
	"out:s"=>\$Out,
	"min:i"=>\$Min,
	"mio:i"=>\$Mio,
	"cpt:i"=>\$Cpt,
	"help"=>\$Help
);

die `pod2text $0` if (!defined $For or !defined $Rev or !defined $Out or defined $Help);

$Cpt ||= 6;
$Min ||= 151;
$Mio ||= 35;

`$Bin/pear -f $For -r $Rev -o $Out -v $Mio -n $Min -s 1 -j $Cpt >pear.log`;

open IN, "$Out\.assembled.fastq" || die $!;
open OUT, ">$Out" || die $!;
my $num=1;
while(<IN>){
	chomp;
	chomp(my $seq=<IN>);
	<IN>;
	<IN>;
	print OUT ">$num\n$seq\n";
	$num++;
}
close OUT;
close IN;
`rm "$Out\.assembled.fastq" "$Out\.discarded.fastq"`;
print "all done\n";
