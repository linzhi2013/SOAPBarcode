#!/usr/bin/env python3
import sys
import os
import re
from Bio import SeqIO


def filter_by_size(in_f=None, size_cutoff=0, out_f=None):
    # >AATGC_265;size=11942
    fh_1 = open(out_f, 'w')
    fh_2 = open(out_f + ".lowSize", 'w')
    for rec in SeqIO.parse(in_f, 'fasta'):
        m = re.search(r'size\=(\d+)', rec.id)
        size = int(m.group(1))
        if size >= size_cutoff:
            SeqIO.write(rec, fh_1, 'fasta')
        else:
            SeqIO.write(rec, fh_2, 'fasta')


def main():
    usage = '''
python3 {} <in.fas> <size_cutoff> <out.fas>

Filter out the sequences with size < size_cutoff.

'''.format(sys.argv[0])

    if len(sys.argv) != 4:
        sys.exit(usage)

    in_f, size_cutoff, out = sys.argv[1:4]
    size_cutoff = int(size_cutoff)
    filter_by_size(in_f, size_cutoff, out)


if __name__ == '__main__':
    main()