#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# use FindBin qw($Bin $Script);
# use PerlIO::gzip;

=head1 Description

        To remove the sequences with <= 0X-depth (by default) sites.

        Usage: perl barcode_translation.pl <parameter>

        --fas           Input Fasta file
        --depth         Depth cutoff. [0]
        --fq            Fastq file, can be used multiple times.
        --out           "prefix of the output file"
        --bwa           Path to 'bwa' executable. [bwa]
        --samtools      Path to 'samtools' executable. [samtools]
        --thread        Thread number. [2]
        --help          "print out this information"

=cut



my ($in_fas,$depth_cutoff, @fqfiles, $prefix, $bwa, $samtools, $thread, $Help);

GetOptions(
        "fas:s"=>\$in_fas,
        "depth:f"=>\$depth_cutoff,
        "fq:s"=>\@fqfiles,
        "out:s"=>\$prefix,
        "bwa:s"=>\$bwa,
        "samtools:s"=>\$samtools,
        "thread:i"=>\$thread,
        "help"=>\$Help
);

die `pod2text $0` if ($Help || !defined ($in_fas) || !defined ($prefix) );
$depth_cutoff = 0 if (!defined $depth_cutoff);
$bwa = "bwa" if (!defined $bwa);
$samtools = "samtools" if (!defined $samtools);
$thread = 2 if (!defined $thread);

# my $fq_files = join(' ', @fqfiles);

# bwa index

my $cmd_index = "$bwa index $in_fas";
print "$cmd_index\n";
# system("$cmd_index") == 0 or die "Command Failed:\n$cmd_index\n";

# bwa sampe
my $mapped_bam = "$prefix.mapped.bam";
my $cmd_bwa = "$bwa mem -t $thread $in_fas @fqfiles | $samtools view -h -b -F 4 | $samtools sort -t $thread -o $mapped_bam ";
print "$cmd_bwa\n";
# system("$cmd_bwa") == 0 or die "Command Failed:\n$cmd_bwa\n";

# get depth
my $depth_file = "$prefix.depth";
# tab-separated
# reference name, position, and coverage depth.
my $cmd_depth = "$samtools depth -a -a -o $depth_file $mapped_bam";
print "$cmd_depth\n";
# system("$cmd_depth") == 0 or die "Command Failed:\n$cmd_depth\n";


# filter sequences by depth
open DE, "<$depth_file" or die $depth_file;
my %hash;
while (<DE>) {
    chomp;
    my ($seqid, $position, $coverage) = split /\t/;
    if ($coverage <= $depth_cutoff){
        $hash{$seqid} = 0;
    }
}
close DE;

my @discard_seqids = keys %hash;
print "Discard these sequences:\n@discard_seqids\n";


open IN, "<$in_fas" or die $in_fas;
my $outfile = "$prefix.depth-lt-${depth_cutoff}X";
open OUT, ">$outfile" or die $outfile;

$/="\>"; <IN>; $/="\n";
while (<IN>) {
    chomp;
    $/="\>";
    chomp(my $seq = <IN>);
    $/="\n";
    $seq=~s/\n//g;

    my $seqid = (split /\s+/)[0];
    unless (exists $hash{$seqid}) {
        print OUT ">$_\n$seq\n";
    }
}

close IN;
close OUT;

