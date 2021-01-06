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
        --view_opt      options for `samtools view` during SAM filtering. ["-h -b -f 0x2"]
        --thread        Thread number. [2]
        --help          "print out this information"

=cut



my ($in_fas,$depth_cutoff, @fqfiles, $prefix, $bwa, $samtools, $samtools_view_opt, $thread, $Help);

GetOptions(
        "fas:s"=>\$in_fas,
        "depth:i"=>\$depth_cutoff,
        "fq:s"=>\@fqfiles,
        "out:s"=>\$prefix,
        "bwa:s"=>\$bwa,
        "samtools:s"=>\$samtools,
        "view_opt:s"=>\$samtools_view_opt,
        "thread:i"=>\$thread,
        "help"=>\$Help
);

die `pod2text $0` if ($Help || !defined ($in_fas) || !defined ($prefix) );
$depth_cutoff = 0 if (!defined $depth_cutoff);
$bwa = "bwa" if (!defined $bwa);
$samtools = "samtools" if (!defined $samtools);
$samtools_view_opt = "-h -b -f 0x2" if (!defined $samtools_view_opt);
$thread = 2 if (!defined $thread);

# my $fq_files = join(' ', @fqfiles);

# bwa index

my $cmd_index = "$bwa index $in_fas";
print "$cmd_index\n";
system("$cmd_index") == 0 or die "Command Failed:\n$cmd_index\n";

# bwa sampe
my $mapped_bam = "$prefix.mapped.bam";
my $cmd_bwa = "$bwa mem -t $thread $in_fas @fqfiles | $samtools view $samtools_view_opt | $samtools view -h -b -F 3840 | $samtools sort -t $thread -o $mapped_bam ";
print "$cmd_bwa\n";
system("$cmd_bwa") == 0 or die "Command Failed:\n$cmd_bwa\n";


# Flag 3840 meaning:
# not primary alignment (0x100)
# read fails platform/vendor quality checks (0x200)
# read is PCR or optical duplicate (0x400)
# supplementary alignment (0x800) # see https://yulijia.net/en/bioinformatics/2015/12/21/Linear-Chimeric-Supplementary-Primary-and-Secondary-Alignments.html#fn:5
#
# *Warning: Flag(s) and 0x8 cannot be set when read is not paired
#                         
# get depth
my $depth_file = "$prefix.depth";
# tab-separated
# reference name, position, and coverage depth.
my $cmd_depth = "$samtools depth -q 20 -Q 30 -s -a -a -o $depth_file $mapped_bam";
print "$cmd_depth\n";
system("$cmd_depth") == 0 or die "Command Failed:\n$cmd_depth\n";


# delete index files
my $cmd_del_index = "rm -rf $in_fas.bwt $in_fas.pac $in_fas.ann $in_fas.amb $in_fas.sa";
print "$cmd_del_index\n";
system("$cmd_del_index") == 0 or die "Command Failed:\n$cmd_del_index\n";

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

