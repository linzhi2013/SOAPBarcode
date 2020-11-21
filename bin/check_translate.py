#!/usr/bin/env python3
import sys
import os
import re
from Bio import SeqIO


def filter_seq_by_translation(rec=None, codon_table=2, frames):
    wi


def main():
    usage = '''Filter out the sequences with stop codons. 

    python3 {} <in.fa> <out.fa> <codon_table> <frame1> [frame2] [frameX] ...'''.format(sys.argv[0])
    if len(sys.argv) < 4:
        sys.exit(usage)

    in_fa, out_fa, codon_table, *frames = sys.argv[1:]

    

if __name__ == '__main__':
    main()
