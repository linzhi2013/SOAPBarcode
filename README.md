# SOAPBarcode

Now move to Github from https://sourceforge.net/projects/metabarcoding/.

# Related publications


Hao M, Jin Q, Meng G, Yang C, Yang S, Shi Z, Tang M, Liu S, Li Y, Li J, Zhang D. Using full-length metabarcoding and DNA barcoding to infer community assembly for speciose taxonomic groups: a case study. Evolutionary Ecology. 2020 Sep 5:1-26. DOI: https://doi.org/10.1007/s10682-020-10072-y

    
# History versions

Can be found on https://sourceforge.net/projects/metabarcoding/.

# Citation

Liu S, Li Y, Lu J, Su X, Tang M, Zhang R, Zhou L, Zhou C, Yang Q, Ji Y, Yu DW. SOAPBarcode: revealing arthropod biodiversity through assembly of Illumina shotgun sequences of PCR amplicons. Methods in Ecology and Evolution. 2013 Dec;4(12):1142-50. DOI: https://doi.org/10.1111/2041-210X.12120

# Changelog

**version 4.4**

modified by Guanliang Meng: 1) upgrade BWA '0.5.9-r16' to
    '0.7.17-r1198-dirty'; 2) filter contigs with internal stop codons
    (translate 1,2,3 frames by default); 3) filter contigs with specific
    sequencing depth by mapping SLS against contigs. **Warning: The two new
    features DO NOT work for samples without barcode sequences ahead of the
    primers.**

**version 4.3**

modified by Guanliang Meng: 

1. use Parallel::ForkManager to
    deal with multi-samples at the same time;

2. support running tasks on localhost and SGE cluster;

3. support resumption of running tasks.

**version 4.2**

modified by Guanliang Meng:

1. use Log::Log4perl and
    Log::Dispatch for logging.

**version 4.1**

modified by Guanliang Meng:

1. add citation infor.

**version 4.0**

modified by Guanliang Meng:

1. change usearch to vsearch;

2. small bug fixed (partially) for reading a lib file (-lib option).

**version 3.0**

modified by Shanlin Liu:

1. protein express check can be skipped without
    setting -int and -pro n;

2. multiple primers with index ahead can be
    proccessed batching, however, if index exists at only one end of primers
    (e.g. former), it should be tackled one by one with parameter -len
    setted as 0;

3. PEAR was introduced as an option for shotgun reads
    overlap

