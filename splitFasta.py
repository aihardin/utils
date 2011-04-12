#! /usr/bin/env python
import os
import sys
from Bio import SeqIO

def main():
    if len(sys.argv) < 2:
        print 'usage: script.py fastq_file'
        print 'fastq file outputed by samtools.pl does not keep everything on the same line'
        sys.exit()
        
    fasta_file = sys.argv[1]
    chrs = ['2R', '2L', '3R', '3L', 'X']
    seqs = SeqIO.parse(open(fasta_file), format='fasta')
    for seq in seqs:
        if seq.id not in chrs: continue
        name = fasta_file.split('.')[0]
        name += '_'+seq.id+'.fasta'
        print name
        SeqIO.write(seq, open(name, 'w'), format='fasta')

if __name__ == '__main__':main()
