#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use FindBin qw($Bin $Script);
# use PerlIO::gzip;

=head1 Description

        To remove the sequences of all frames with stop codons or unknown codons.

        Usage: perl barcode_translation.pl <parameter>

        --fas           Input Fasta file
        --frame         Frames to translate. [1,2,3]
        --code          Genetic code [5]
        --out           "prefix of the output file"
        --seqkit        Path to 'seqkit' executable. [seqkit]
        --help          "print out this information"

=cut



my ($in_fas,$frame,$genetic_code,$prefix, $seqkit, $Help);

GetOptions(
        "fas:s"=>\$in_fas,
        "frame:s"=>\$frame,
        "code:i"=>\$genetic_code,
        "out:s"=>\$prefix,
        "seqkit:s"=>\$seqkit,
        "help"=>\$Help
);

die `pod2text $0` if ($Help || !defined ($in_fas) || !defined ($prefix) );
$frame = "1,2,3" if (!defined $frame);
$genetic_code = 5 if (!defined $genetic_code);
$seqkit = "seqkit" if (!defined $seqkit);

my $frame_count = split /,/, $frame;

# seqkit translate --allow-unknown-codon --clean --frame 1,2,3 --transl-table 5 --out-file a HeBi.cds

my $protein_file = "$prefix.translated";
my $cmd1 = "$seqkit translate --allow-unknown-codon --clean --frame $frame --transl-table $genetic_code --out-file $protein_file $in_fas";

print "$cmd1\n";
system("$cmd1") == 0 or die "Command Failed:\n$cmd1\n";

my %protein_hash;
open (PRO, "<$protein_file")  or die $protein_file;
$/="\>"; <PRO>; $/="\n";
while (<PRO>) {
    chomp;
    $/="\>";
    chomp(my $protein = <PRO>);
    $/="\n";
    $protein=~s/\n//g;

    my $ti = $_;
    $ti =~ s/\_frame\=\d+.*$//g;
    unless (exists $protein_hash{$ti}) {
        $protein_hash{$ti} = 0;
    }

    if ($protein=~/X|\*/) {
        $protein_hash{$ti} += 1;
    }

}
close PRO;


my %cds_hash;
open (CDS, "<$in_fas")  or die $in_fas;
$/="\>"; <CDS>; $/="\n";
while (<CDS>) {
    chomp;
    $/="\>";
    chomp(my $seq = <CDS>);
    $/="\n";
    $seq=~s/\n//g;

    my $ti = $_;
    $ti =~ s/^(\S+).*$/$1/g;
    unless (exists $cds_hash{$ti}) {
        $cds_hash{$ti} = ">$_\n$seq";
    }

}
close CDS;


my $clean_cds_file = "$prefix.clean.fas";
my $error_cds_file = "$prefix.X.fas";
open (OUT, ">$clean_cds_file") or die $clean_cds_file;
open (OUT2, ">$error_cds_file") or die $error_cds_file;

for my $ti (keys %protein_hash) {
    my $stop_count = $protein_hash{$ti};
    my $cds_seq = $cds_hash{$ti};

    if ($stop_count < $frame_count) {
        print OUT "$cds_seq\n";
    }else{
        print OUT2 "$cds_seq\n";
    }

}
close OUT;
close OUT2;



