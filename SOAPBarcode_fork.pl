#! usr/bin/perl -w

=head1 Desription

This program is attempt to assemble standard COI barcode region with two libraries, of which one is full length, another is shotgun.

=head1 Version

version 4.3
modified by Guanliang Meng: 1) use Parallel::ForkManager to deal with
multi-samples at the same time;
2) support local and SGE for each sample assembly.

version 4.2
modified by Guanliang Meng: 1) use Log::Log4perl and Log::Dispatch for logging.

version 4.1
modified by Guanliang Meng: 1) add citation infor.

version 4.0
modified by Guanliang Meng: 1) change usearch to vsearch;
          2) small bug fixed (partially) for reading a lib file (-lib option).

version 3.0
modified: 1) protein express check can be skipped without setting -int and -pro n;
      2) multiple primers with index ahead can be proccessed batching, however, if index
         exists at only one end of primers (e.g. former), it should be tackled one by one
         with parameter -len setted as 0;
      3) PEAR was introduced as an option for shotgun reads overlap

Scripted by Jiangliang Lu && Shanlin LIU, any question can contact via liushanlin@genomics.cn

=head1 Usage

perl SOAPBarcode.pl <parameter>


---------------| the sequence information of two libs |-----------------

    -lib    lib file include information of all the raw data. Empty lines after
        each "[LIB]" for each library will lead to bugs!!!

-------------------------| denoise parameter|--------------------------

    -pri    the primer set include barcode sequnces ahead
    -len    the length of the barcode region, set as 0 when no index (0)
    -qac    the maximun B quality allowed for each seq(5)
    -abc    the minimum abandance of the uniq reads(2)
    -ucf    the similarity cutoff when clustering denoise(0.98)
    -int    the hydrophily interval of each site generated by Pro_C,
                formated in "site max min"
    -pro    protein coding gene expression check, (y|n)
    -mpr    the maximum number of mismatch check for primers including indexs ahead (1),
    -oop    choose the overlap programs, 1 for COPE; 2 for PEAR (1),
    -osc    the similarity cutoff for shotgun reads overlaping for COPE (0.95)

-----------------------| assembly parameter |---------------------------

    -lmk    the maximum length of the kmer (115)
    -lsk    the minimum length of the kmer (95)
    -kin    the kmer interval (10)
    -clk    the lower frequency kmer cutoff (0.1)
    -clb    the lower support branch cutoff (0.1)
    -lms    the maximun extension length (660)
    -lss    the minimun extension length (450)
    -cpt    the CPU number allowed for each assembly task (8)

--------------------------| other parameter |---------------------------

    -ucs    the similarity cutoff of the OTU generation (0.98)
    -out    prefix of name of output file
    -ptask  the maximum paralle tasks allowed (1)
    -sge    submit the assembly tasks to SGE cluster to run, default local.
    -qsubt  qsub command template, default 'qsub -cwd -l vf={vf} -pe smp {cpu}'.
            You always have the variables '{vf}' and '{cpu}'.
    -avf    the '{vf}' value for each assembly task ('50G')
    -resume Resume the run.
    -tmpdir tmp directory path (./tmp.soapbarcode)
    -help   print out this information
    -debug  Debug mode

-----------------------------| Example |-------------------------------

commond:
perl SOAPBarcode.pl -lib test.lib -pri primer.fasta -int interval -pro y -out test -oop 2 -len 5 -mpr 0

=head1 Citation

    Liu S, Li Y, Lu J, Su X, Tang M, Zhang R, Zhou L, Zhou C, Yang Q, Ji Y, Yu DW. 
    SOAPBarcode: revealing arthropod biodiversity through assembly of Illumina
    shotgun sequences of PCR amplicons. Methods in Ecology and Evolution. 2013 Dec;4(12):1142-50.
    DOI: https://doi.org/10.1111/2041-210X.12120

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use Log::Log4perl qw(get_logger :levels);
use Log::Dispatch;
use Parallel::ForkManager;
use Env qw(PATH PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT);
# use Cwd;
use Cwd qw(getcwd abs_path);

my $workdir = getcwd;

my $logger = get_logger("SOAPBarcode");

my ($Debug, $Parallel_task, $submit_sge, $resume, $avf, $qsubt, $tmpdir);
my ($Lib,$Ffq,$Bfq,$Fsfq,$Bsfq,$Primer,$Bcut,$Interval,$Out,$Help,$Len,$Minimum,$Mpr,$Pro,$Oop);
my ($Osc,$Lmk,$Lsk,$Kin,$Clk,$Clb,$Lms,$Lss,$Cpt,$Ucs,$Ucf);
GetOptions(
    "lib:s"=>\$Lib,
    "pro:s"=>\$Pro,
    "pri:s"=>\$Primer,
    "qac:i"=>\$Bcut,
    "len:i"=>\$Len,
    "abc:i"=>\$Minimum,
    "int:s"=>\$Interval,
    "mpr:i"=>\$Mpr,
    "osc:i"=>\$Osc,
    "oop:i"=>\$Oop,
    "lmk:i"=>\$Lmk,
    "lsk:i"=>\$Lsk,
    "kin:i"=>\$Kin,
    "clk:i"=>\$Clk,
    "clb:i"=>\$Clb,
    "lms:i"=>\$Lms,
    "lss:i"=>\$Lss,
    "cpt:i"=>\$Cpt,
    "ucf:i"=>\$Ucf,
    "ucs:i"=>\$Ucs,
    "out:s"=>\$Out,
    "debug"=>\$Debug,
    "ptask:i"=>\$Parallel_task,
    "sge"=>\$submit_sge,
    "qsubt"=>\$qsubt,
    "avf"=>\$avf,
    "resume"=>\$resume,
    "tmpdir:s"=>\$tmpdir,
    "help"=>\$Help
);

if(defined $Debug) {
    $logger->level($DEBUG);
}else{
    $logger->level($INFO);
}

my $appender = Log::Log4perl::Appender->new(
    "Log::Dispatch::File",
    filename => "SOAPBarcode.log",
    mode     => "append",
 );
 $logger->add_appender($appender);
# DEBUG
# INFO
# WARN
# ERROR
# FATAL


$Parallel_task ||= 1;
if ($Parallel_task == 1){
    $Parallel_task = 0;
}

if ($Parallel_task < 0) {
    $logger->error("-ptask must be >= 0");
}

$qsubt = 'qsub -cwd -l vf={vf} -pe smp {cpu}' if (!defined $qsubt);
$avf ||= "50G";
$tmpdir = './tmp.soapbarcode' if (!defined $tmpdir);

# for testing
if (-e "$tmpdir") {
    $logger->info("$tmpdir exists\n!")
}else{
    system("mkdir -p $tmpdir",)
}
$tmpdir = abs_path($tmpdir);

# run_cmd($logger,
#    "Creating $tmpdir directory",
#    "mkdir -p $tmpdir",
#    '');

$Bcut=5 if (!defined $Bcut);
$Len ||= 0;
$Minimum ||=2;
$Mpr=1 if (!defined $Mpr);
$Osc ||=0.95;
$Oop ||=1;
$Lmk ||=115;
$Lsk ||=91;
$Kin ||=10;
$Clk ||=0.1;
$Clb ||=0.1;
$Lms ||=660;
$Lss ||=450;
$Cpt ||=8;
$Ucs ||=0.98;
$Ucf ||=0.98;
die `pod2text $0` if ($Help || !defined $Lib || (defined $Interval && ($Pro eq "n"))|| !defined $Primer || !defined $Out);

if (!$Pro eq "y" and !$Pro eq "n"){
    $logger->error("-pro parameter should be y or n for protein expression check option");
}
if (!$Oop == 1 and !$Oop == 2){
    $logger->error("-oop parameter should be setted as 1 for COPE and as 2 for PEAR");
}

open LIB, $Lib || die $!;
my (@in_si, @lrl, @fastq, $in_sif,$in_sis,$lrlf,$lrls);
while(<LIB>){
    chomp;
    next if(/^\s*$/); # mgl. fixed a bug if the lib file has emptt lines.
    my ($a, $b, $c, $d);
    if (/\[LIB\]/){
        chomp($a=<LIB>);
        chomp($b=<LIB>);
        chomp($c=<LIB>);
        chomp($d=<LIB>);
        $logger->error("wrong lib file, makesure the second line represents the insertsize") unless ($a=~s/insertsize=//);
        $logger->error("wrong lib file, makesure the second line represents the readlength") unless ($b=~s/readlength=//);
        $logger->error("wrong lib file, makesure the second line represents the fastq 1") unless ($c=~s/fq1=//);
        $logger->error("wrong lib file, makesure the second line represents the fastq 2") unless ($d=~s/fq2=//);
        push @in_si, $a;
        push @lrl, $b;
        push @fastq, $c;
        push @fastq, $d;
    }
}
$logger->error("wrong lib file, make sure the lib file include two libs with pair end reads") unless (@in_si==2 && @lrl==2 && @fastq==4);
if ($in_si[0] >= $in_si[1]){
    $in_sif=$in_si[0];
    $in_sis=$in_si[1];
    $lrlf=$lrl[0];
    $lrls=$lrl[1];
    $Ffq=$fastq[0];
    $Bfq=$fastq[1];
    $Fsfq=$fastq[2];
    $Bsfq=$fastq[3];
}else{
    $in_sif=$in_si[1];
    $in_sis=$in_si[0];
    $lrlf=$lrl[1];
    $lrls=$lrl[0];
    $Ffq=$fastq[2];
    $Bfq=$fastq[3];
    $Fsfq=$fastq[0];
    $Bsfq=$fastq[1];
}

########################get full length reads set###########################
my $f_out="$Out"."_f";
my @flist;
my $binpath= $Bin;
if ($Mpr == 0){
    run_task($logger,
        "Assign_FLS2DiffSamples_By_Tag",
        "Assign FLS data to different samples:",
        "perl $binpath/bin/extract.perfect.pl -fq1 $Ffq -fq2 $Bfq -pri $Primer -out $f_out -len $Len -mis $Mpr -qua $Bcut >>log",
        "$Ffq + $Bfq --> $f_out\_list",
        $resume,
        $submit_sge,
        $qsubt,
        "500M",
        1,
        $workdir,
        $tmpdir);
}else{
    run_task($logger,
        "Assign_FLS2DiffSamples_By_Tag",
        "Assign FLS data to different samples:",
        "perl $binpath/bin/extract.v1.pl -fq1 $Ffq -fq2 $Bfq -pri $Primer -out $f_out -len $Len -mis $Mpr -qua $Bcut >>log",
        "$Ffq + $Bfq --> $f_out\_list",
        $resume,
        $submit_sge,
        $qsubt,
        "500M",
        1,
        $workdir,
        $tmpdir);
}
$logger->info("fq1=$Ffq\nfq2=$Bfq\nreads has been assigned\n", "FLS files for all samples are in $f_out\_list\n");

my $pm_fls_denoise = Parallel::ForkManager->new($Parallel_task);
# Setup a callback for when a child finishes up so we can
# get it's exit code
$pm_fls_denoise->run_on_finish( sub {
    my ($pid, $exit_code, $ident) = @_;
    print "** $ident just got out of the pool ".
      "with PID $pid and exit code: $exit_code\n";
});

$pm_fls_denoise->run_on_start( sub {
    my ($pid, $ident)=@_;
    print "** $ident started, pid: $pid\n";
});

$pm_fls_denoise->run_on_wait( sub {
    print "** Have to wait for one children ...\n"
  },
  0.5
);

open FLI, "$f_out\_list" || die $!;
#open FEN, ">$fendlist" || die $!;
my $full_out_pct=(100-$Ucf*100);

FLS_DENOISE_LOOP:
while(<FLI>){
    chomp;
    # Forks and returns the pid for the child:
    my $pid = $pm_fls_denoise->start($_) and next FLS_DENOISE_LOOP;

    my $dupout="$_".".dup";
    run_task($logger,
        "FLS_Dereplication.$dupout",
        "Dereplication of FLS data:",
        "$binpath/bin/vsearch --derep_fulllength $_ --output $dupout --sizeout --strand both --minuniquesize 2",
        "$_ --> $dupout",
        $resume,
        $submit_sge,
        $qsubt,
        "500M",
        1,
        $workdir,
        $tmpdir);

    if ($Pro eq "y") {
        my $proout="$_".".pro";
        run_task($logger,
            "FLS_PCG_expression_check.$proout",
            "Protein coding gene expression check for FLS data:",
            "perl $binpath/bin/Pro_C.pl -lib f -int $Interval -fas $dupout -out $proout",
            "$dupout --> $proout",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);
        my $otuout="$_".".otu";
        run_task($logger,
            "FLS_CLUSTER.$otuout",
            "Clustering of FLS data:",
            "$binpath/bin/vsearch --cluster_size $proout --threads 1 --sizein --sizeout --id $Ucf --centroids $otuout",
            "$proout --> $otuout",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);

        my $endout="$_".".end";
        $logger->info("Split the connected FLS reads to read1 and read2:\n", "$otuout --> $endout");
        open TEI, "$otuout" || die $!;
        open TEO, ">$endout" || die $!;
        $/="\>";<TEI>;$/="\n";
        while(my $ti=<TEI>){
            chomp($ti);
            $/="\>";
            chomp(my $seq = <TEI>);
            $seq=~s/\n//g;
            my @a=split /NNN/,$seq;
            $logger->error("the delimiter is wrong, reads can't be seperated into two part.\nRead name: $ti\nFilename: $otuout") unless (@a == 2);
            my $endlenf=length $a[0];
            my $endlenr=length $a[1];
            if ($endlenf > $Lmk and $endlenr> $Lmk){
                print TEO ">$ti\n$a[0]\n$a[1]\n";
            }
            $/="\n";
        }
        close TEI;
        close TEO;
#       print FEN "$endout\n";
        # push @flist,$endout; # mgl: no use in paralle task
        run_cmd($logger,
            "Removing intermediate files: $dupout $proout $otuout",
            "rm $dupout $proout $otuout", "");
    }elsif ($Pro eq "n") {
        my $otuout="$_".".otu";
        run_task($logger,
            "FLS_CLUSTER.$otuout",
            "Skip protein coding gene expression check for FLS data...\nClustering of FLS data:",
            "$binpath/bin/vsearch --cluster_size $dupout --threads 1 --sizein --sizeout --id $Ucf --centroids $otuout",
            "$dupout --> $otuout",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);

        my $endout="$_".".end";
        $logger->info("Split the connected FLS reads to read1 and read2:\n", "$otuout --> $endout");
        open TEI, "$otuout" || die $!;
        open TEO, ">$endout" || die $!;
        $/="\>";<TEI>;$/="\n";
        while(my $ti=<TEI>){
            chomp($ti);
            $/="\>";
            chomp(my $seq = <TEI>);
            $seq=~s/\n//g;
            my @a=split /NNN/,$seq;
            $logger->error("the delimiter is wrong, reads can't be seperated into two part.\nRead name: $ti\nFilename: $otuout") unless (@a == 2);

            my $endlenf=length $a[0];
            my $endlenr=length $a[1];
            if ($endlenf > $Lmk and $endlenr> $Lmk){
                print TEO ">$ti\n$a[0]\n$a[1]\n";
            }
            $/="\n";
        }
        close TEI;
        close TEO;
#       print FEN "$endout\n";
        # push @flist,$endout; # mgl: no use in paralle task
        run_cmd($logger,
            "Removing intermediate files: $dupout $otuout",
            "rm $dupout $otuout", "");
    }else{
        $logger->error("-pro parameter should be y or n for protein expression check option");
    }

    $pm_fls_denoise->finish;
}
$logger->info("FLS_DENOISE_LOOP: Waiting for Children...\n");
$pm_fls_denoise->wait_all_children;
$logger->info("FLS_DENOISE_LOOP: Everybody is out of the pool!\n");

close FLI;

@flist = glob("$f_out*.fasta.end"); # mgl: get all *.fasta.end files.

########################get full length reads set###########################

#---------------------------------------------------------------#
#---------------------------------------------------------------#

##########################get shotgun reads set##############################
my $Out4="$Out"."."."shotgun";
my ($cmrl,$cmru);
my $expect=2*$lrls-$in_sis;
$cmru=($lrls-1);
if ($expect>=60){$cmrl=int($expect/3)}else{$cmrl=20};
if ($Oop ==1 ){
    run_task($logger,
        "SLS_connect_R1_R2_CMR",
        "Connect read1 and read2 for SLS data:",
        "$binpath/bin/cmr -a $Fsfq -b $Bsfq -o $Out4 -2 $Out.1.left -3 $Out.2.left -l $cmrl -u $cmru -c $Osc >overlap.log",
        "$Fsfq + $Bsfq --> $Out4",
        $resume,
        $submit_sge,
        $qsubt,
        "500M",
        1,
        $workdir,
        $tmpdir);
}elsif ($Oop ==2){
    my $pear_min=$lrls+1;
    my $pear_mio=$cmrl;
    run_task($logger,
        "SLS_connect_R1_R2_PEAR",
        "Connect read1 and read2 for SLS data:",
        "perl $binpath/bin/pear_overlap.pl -for $Fsfq -rev $Bsfq -out $Out4 -min $pear_min -mio $pear_mio -cpt $Cpt",
        "$Fsfq + $Bsfq --> $Out4",
        $resume,
        $submit_sge,
        $qsubt,
        "500M",
        1,
        $workdir,
        $tmpdir);

}else{
    $logger->error("the overlap program can be setted as 1 for COPE and as 2 for PEAR");
}
$logger->info("$Fsfq and $Bsfq: the shotgun reads have been overlapped, result file: $Out4");

my $s_out="$Out"."_s";
if ($Mpr ==0){
    run_task($logger,
        "Assign_SLS2DiffSamples_By_Tag",
        "Assign SLS data to different samples:",
        "perl $binpath/bin/shotgun_assign.perfect.pl -pri $Primer -fas $Out4 -out $s_out -mis $Mpr -len $Len",
        "$Out4 --> $s_out\_list",
        $resume,
        $submit_sge,
        $qsubt,
        "500M",
        1,
        $workdir, 
        $tmpdir);
}else{
    run_task($logger,
        "Assign_SLS2DiffSamples_By_Tag",
        "Assign SLS data to different samples:",
        "perl $binpath/bin/shotgun_assign.pl -pri $Primer -fas $Out4 -out $s_out -mis $Mpr -len $Len",
        "$Out4 --> $s_out\_list",
        $resume,
        $submit_sge,
        $qsubt,
        "500M",
        1,
        $workdir,
        $tmpdir);
}
$logger->info("Shotgun reads ($Out4) has been assigned\n", "SLS files for all samples are in $s_out\_list\n");


my $pm_sls_denoise = Parallel::ForkManager->new($Parallel_task);
# Setup a callback for when a child finishes up so we can
# get it's exit code
$pm_sls_denoise->run_on_finish( sub {
    my ($pid, $exit_code, $ident) = @_;
    print "** $ident just got out of the pool ".
      "with PID $pid and exit code: $exit_code\n";
});

$pm_sls_denoise->run_on_start( sub {
    my ($pid, $ident)=@_;
    print "** $ident started, pid: $pid\n";
});

$pm_sls_denoise->run_on_wait( sub {
    print "** Have to wait for one children ...\n"
  },
  0.5
);

open MLI, "$s_out\_list" || die $!;
my @slist;
my $s_bwa;
#open MDU, ">$mdout" || die $!;
SLS_DENOISE_LOOP:
while (<MLI>) {
    chomp;

    # Forks and returns the pid for the child:
    my $pid = $pm_sls_denoise->start($_) and next SLS_DENOISE_LOOP;

    my $sortout="$_".".sort";
    if ($Pro eq "y") {
        my $proout="$_".".pro";
        run_task($logger,
            "SLS_PCG_expression_check.$sortout",
            "Protein coding gene expression check for SLS data:",
            "perl $binpath/bin/Pro_C.pl -lib s -fas $_ -out $proout",
            "$_ --> $proout",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);
        run_cmd($logger,
            "Renaming $proout to $sortout",
            "mv $proout $sortout", "");
#       `$binpath/bin/usearch -sortbylength $proout -output $sortout`;
#       `rm $proout`;
    }elsif ($Pro eq "n") {
        $logger->info("Skip protein coding gene expression check for SLS data\n", "Renaming $_ to $sortout");
        `mv $_ $sortout`;
#       `$binpath/bin/usearch -sortbylength $_ -output $sortout`;
    }else{
        $logger->error("-pro parameter should be y or n for protein expression check option");
    }

    my $dupout="$_".".dup";
    run_task($logger,
        "FLS_Dereplication.$dupout",
        "Dereplication of SLS data:",
        "perl $binpath/bin/rmdupctg.pl $sortout $dupout",
        "$sortout --> $dupout",
        $resume,
        $submit_sge,
        $qsubt,
        "500M",
        1,
        $workdir,
        $tmpdir);

#   print MDU "$dupout\n";
    # push @slist, $dupout; # mgl: no use in paralle task
    if ($Len==0){
        $s_bwa = $sortout
    }else{
        $logger->info("Removing $sortout");
        `rm $sortout`
    }

    $pm_sls_denoise->finish;
}
$logger->info("SLS_DENOISE_LOOP: Waiting for Children...\n");
$pm_sls_denoise->wait_all_children;
$logger->info("SLS_DENOISE_LOOP: Everybody is out of the pool!\n");

@slist = glob("$s_out*.fasta.dup"); # mgl: get all *.fasta.dup files.

close MLI;
$logger->info("The shotgun reads have been denoised, result files: @slist");

$logger->info("Removing $s_out\_list $f_out\_list");
 `rm "$s_out\_list" "$f_out\_list"`;
##########################get shotgun reads set##############################

#---------------------------------------------------------------#
#---------------------------------------------------------------#

##########################      asembly        ##############################

my %component;
for my $i (0..$#flist){
    my @a=split /\_|\./,$flist[$i];
    $component{$a[2]}[0]=$flist[$i]; # mgl: 0: *.end files from FLS
}
for my $i (0..$#slist){
    my @a=split /\_|\./,$slist[$i];
    if (exists $component{$a[2]}){
        push @{$component{$a[2]}},$slist[$i]; # mgl: 1: *.dup files from SLS
    }
}

my $pm_assemble = Parallel::ForkManager->new($Parallel_task);
# Setup a callback for when a child finishes up so we can
# get it's exit code
$pm_assemble->run_on_finish( sub {
    my ($pid, $exit_code, $ident) = @_;
    print "** $ident just got out of the pool ".
      "with PID $pid and exit code: $exit_code\n";
});

$pm_assemble->run_on_start( sub {
    my ($pid, $ident)=@_;
    print "** $ident started, pid: $pid\n";
});

$pm_assemble->run_on_wait( sub {
    print "** Have to wait for one children ...\n"
  },
  300
);

BARCODE_ASSEMBLE:
for my $key (keys %component) {
    if (@{$component{$key}}!=2) {
        $logger->warn("index $key only generate one lib, skip this!");
        next;
    }

    # Forks and returns the pid for the child:
    my $pid = $pm_assemble->start($key) and next BARCODE_ASSEMBLE;

    my $libout="$tmpdir/$key\.lib";
    $logger->info("Preparing $libout for sample $key");
    open LIB, ">$libout" || die $!;
    print LIB ">\nf=$workdir/$component{$key}[1]";
    close LIB;
    my $assout="$Out\_$key";

    run_task($logger,
        "barcode_assembly.$key",
        "Let's assemble $key:",
        "$binpath/bin/barcode -e $component{$key}[0] -r $libout -l $Lsk -k $Lmk -o $assout -v $Kin -c $Clk -s $Clb -n $Lss -x $Lms -t $Cpt 1>$tmpdir/$key.barcode.log 2>$tmpdir/$key.barcode.err",
        "Result file: $assout.contig",
        $resume,
        $submit_sge,
        $qsubt,
        $avf,
        $Cpt,
        $workdir,
        $tmpdir);

    $logger->info("Choose the first isoform for each locus:\n $assout.contig --> $assout.contig.F");
    open ASS, "$assout.contig" || die $!;
    open ASF, ">$assout.contig.F" || die $!;
    my %aha;
    $/="\>"; <ASS>; $/="\n";
    while(<ASS>){
        chomp;
        $/="\>";
        chomp(my $seq = <ASS>);
        $/="\n";
        $seq=~s/\n//g;
        my $len=length $seq;
        my @a=split /\t/;
        # >115_1_1        533 # mgl: only choose the first one
        # >115_1_2        815
        # >115_2_1        533
        # >115_3_1        533
        my @b=split /\_/, $a[0];
        $aha{$b[1]}=0 unless(exists $aha{$b[1]});
        if ($b[2]==1){
                $aha{$b[1]}++;
                print ASF ">$_\n$seq\n";
        }elsif($aha{$b[1]}==0){
                $aha{$b[1]}++;
                print ASF ">$_\n$seq\n";
        }else{
                next;
        }
    }
    close ASS;
    close ASF;
#   push @assembled, "$assout.contig.F";
    # below: mgl: no use in paralle task
    # push @{$component{$key}},"$assout.contig.F"; # mgl: 2: *.contig.F files from barcode program

    $pm_assemble->finish;
}

$logger->info("BARCODE_ASSEMBLE: Waiting for Children...\n");
$pm_assemble->wait_all_children;
$logger->info("BARCODE_ASSEMBLE: Everybody is out of the pool!\n");

for my $key (keys %component) {
    my $assout="$Out\_$key";
    push @{$component{$key}},"$assout.contig.F";
}

##########################      asembly        ##############################

#---------------------------------------------------------------#
#---------------------------------------------------------------#

#######   abundance information retrieve and final cluster    #############

$logger->info("Assembly finished!\nLet's obtain the abundance information!");

$full_out_pct=(100-$Ucs*100);
my $com_num=keys %component;

my $pm_abundance = Parallel::ForkManager->new($Parallel_task);
# Setup a callback for when a child finishes up so we can
# get it's exit code
$pm_abundance->run_on_finish( sub {
    my ($pid, $exit_code, $ident) = @_;
    print "** $ident just got out of the pool ".
      "with PID $pid and exit code: $exit_code\n";
});

$pm_abundance->run_on_start( sub {
    my ($pid, $ident)=@_;
    print "** $ident started, pid: $pid\n";
});

$pm_abundance->run_on_wait( sub {
    print "** Have to wait for one children ...\n"
  },
  0.5
);

GET_ABUNDANCE:
for my $key (keys %component){
    next unless (defined $component{$key}[2]);

    # Forks and returns the pid for the child:
    my $pid = $pm_abundance->start($key) and next GET_ABUNDANCE;

    if ($com_num==1 && $Len==0){
        $logger->info("To obtain the abundance information by mapping SLS reads against assembled sequences");

        run_task($logger,
            "$component{$key}[2].bwa.index",
            "Indexing the assembly file: $component{$key}[2]:",
            "$binpath/bin/bwa index $component{$key}[2]",
            "",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);

        run_task($logger,
            "$component{$key}[2].bwa.aln",
            "Maping SLS data $s_bwa to $component{$key}[2]:",
            "$binpath/bin/bwa aln -n 0 -t $Cpt $component{$key}[2] $s_bwa >bwa.sai",
            "",
            $resume,
            $submit_sge,
            $qsubt,
            "2G",
            $Cpt,
            $workdir,
            $tmpdir);

        run_task($logger,
            "$component{$key}[2].bwa.samse",
            "Generating SAM file:",
            "$binpath/bin/bwa samse $component{$key}[2] bwa.sai $s_bwa >bwa.sam",
            "",
            $resume,
            $submit_sge,
            $qsubt,
            "1G",
            1,
            $workdir,
            $tmpdir);

        my $abun_out="$component{$key}[2]"."A";
        run_task($logger,
            "$component{$key}[2].cal.abundance",
            "Calculate average sequencing coverage based on SAM file:",
            "perl $binpath/bin/sam_abundance.pl -fas $component{$key}[2] -sam bwa.sam -out $abun_out",
            "$component{$key}[2] + bwa.sam --> $abun_out",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);

        my $sort_out="$abun_out"."S";
        run_task($logger,
            "$component{$key}[2].barcode.sortbysize",
            "Sort $abun_out by size --> $sort_out",
            "$binpath/bin/vsearch --sortbysize $abun_out --output $sort_out",
            "",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);

        my $cluster_out="$sort_out"."TA";
        run_task($logger,
            "$component{$key}[2].barcode.cluster",
            "Cluster $sort_out into OTUs --> $cluster_out:",
            "$binpath/bin/vsearch --cluster_size $sort_out --threads 1 --sizein --sizeout --id $Ucs --centroids $cluster_out",
            "",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);
        run_cmd($logger,
            "Removing $abun_out $sort_out",
            "rm $abun_out $sort_out", "");

    }elsif($com_num >1 && $Len>0){
        $logger->info("To obtain the abundance information for $component{$key}[2] by extracting 'size' info from $component{$key}[0] of FLS data --> $component{$key}[2]A");
        my %chash;
        open CAS, "$component{$key}[2]" || die $!; # mgl: the *.contig.F files from barcode program
        $/="\>";<CAS>;$/="\n";
        while(<CAS>){
            chomp;
            $/="\>";
            chomp(my $seq=<CAS>);
            $/="\n";
            my @a=split /\_/,$_;
            $chash{$a[1]}[0]=$seq;
        }
        close CAS;
        my $cnum=1;
        open CFL, "$component{$key}[0]" || die $!;
        while(<CFL>){
            next unless (/^>/);
            # >TACCT_179;size=1616
            $chash{$cnum}[1]=$_ if (exists $chash{$cnum});
            $cnum++;
        }
        close CFL;
        my $abun_out="$component{$key}[2]"."A";
        open COU, ">$abun_out" || die $!;
        for my $ckey(keys %chash) {
            print COU "$chash{$ckey}[1]\n$chash{$ckey}[0]";
        }
        close COU;

        my $sort_out="$abun_out"."S";
        run_task($logger,
            "$component{$key}[2].barcode.sortbysize",
            "Sort $abun_out by size --> $sort_out",
            "$binpath/bin/vsearch --sortbysize $abun_out --output $sort_out",
            "",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);

        my $cluster_out="$sort_out"."TA";
        run_task($logger,
            "$component{$key}[2].barcode.cluster",
            "Cluster $sort_out into OTUs --> $cluster_out:",
            "$binpath/bin/vsearch --cluster_size $sort_out --threads 1 --sizein --sizeout --id $Ucs --centroids $cluster_out",
            "",
            $resume,
            $submit_sge,
            $qsubt,
            "500M",
            1,
            $workdir,
            $tmpdir);

#        `rm $abun_out $sort_out`;
    }else{
        $logger->error("there are problems of index length ($Len) and number of assembled reuslts ($com_num)");
    }

    $pm_abundance->finish;
}

$logger->info("GET_ABUNDANCE: Waiting for Children...\n");
$pm_abundance->wait_all_children;
$logger->info("GET_ABUNDANCE: Everybody is out of the pool!\n");


#######   abundance information retrieve and final cluster    #############

$logger->info("all done");

sub run_cmd {
    my ($logger, $prefix_msg, $cmd, $post_msg) = @_;
    $logger->info("$prefix_msg\n$cmd\n$post_msg\n");
    system("$cmd") == 0 or $logger->error("Command Failed:\n$cmd\n");
}


sub run_task{
    my ($logger, $job_iden, $prefix_msg, $cmd, $post_msg, $resume, $submit_sge, $qsubt, $vf, $cpu, $workdir, $tmpdir) = @_;
    my $done_file = "$tmpdir/$job_iden.done";
    my $stderr_file = "$tmpdir/$job_iden.error";
    my $stdout_file = "$tmpdir/$job_iden.log";

    if ($resume and -e $done_file) {
        $logger->info("Use existing result for step $job_iden\n");
        return 0;
    }elsif (-e $done_file) {
        $logger->info("Removing $done_file...");
        system("rm $done_file");
    }

    if ($submit_sge) {
        $logger->info("$prefix_msg\n$cmd\n$post_msg\nThis task will be submitted to SGE\n");
        my $shell_file = "$tmpdir/$job_iden.sh";
        open OUT, ">$shell_file" or die $!;
        # print OUT "#!/usr/bin/bash\n";
        my $message = <<"END_MESSAGE";
############################################
# You may need to adapt this manually!!! ###
############################################
export PATH='$PATH'
export PERL5LIB='$PERL5LIB'
export PERL_LOCAL_LIB_ROOT='$PERL_LOCAL_LIB_ROOT'
export PERL_MB_OPT='$PERL_MB_OPT'
export PERL_MM_OPT='$PERL_MM_OPT'

cd "$workdir"

set -vex
$cmd
touch $done_file
END_MESSAGE

        print OUT $message;
        close OUT;

        $qsubt=~s/\{vf\}/$vf/g;
        $qsubt=~s/\{cpu\}/$cpu/g;
        my $submit_sge = "$qsubt -e $stderr_file -o $stdout_file $shell_file";
        system("$submit_sge") == 0
            or $logger->error("Command Failed:\n$submit_sge\n");

        while (1) {
            sleep 10;
            last if (-e "$done_file");
        }
        $logger->info("Command finished:\n$submit_sge\n");

    }else{
        $logger->info("$prefix_msg\n$cmd\n$post_msg\nThis task will run locally\n");
        system("$cmd") == 0 or $logger->error("Command Failed:\n$cmd\n");
        system("touch $done_file") == 0 or $logger->error("Command Failed:\ntouch $done_file\n");
    }

}


