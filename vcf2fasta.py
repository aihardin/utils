#! /usr/bin/env python
import os
import sys
from Bio import SeqIO

def readFasta(fastafn):
    '''Read a fasta file, return a dictionary keyed by sequence name containing the respective sequences.'''
    fastafh = open(fastafn, 'r')
    
    # parse the FASTA file
    seqs = {}
    for s in list(SeqIO.parse(fastafh, format='fasta')):
        seqs[s.name] = s.seq.tostring()

    return seqs


def parseVCF(snp_file,filter_vcf=True,filter_het=True,write_filtered=None):
    subs = list()
    indels = list()
    #strain = snp_file.split('.')[0]
    if write_filtered:
        filtered_fh = open(write_filtered,'w')
    
    for line in snp_file:
        if line.startswith('#'):
            if write_filtered:
                filtered_fh.write(line)
            if line.startswith('#CHROM'):
                strain = line.split()[9].split('.')[0]
            continue
        if filter_vcf:
            if 'PASS' not in line:
                continue
        if filter_het:
            if not '1/1' in line and not '1|1' in line:
                continue
                
        if write_filtered:
            filtered_fh.write(line)

        p = line.split('\t')
        chr =p[0]
        pos = int(p[1])
        ref = p[3]
        var = p[4]
        genotype = p[9].split(':')[0]
        if len(var) != len(ref):
            indels.append((chr,pos, ref, var, genotype))
        else:
            subs.append((chr, pos, ref, var, genotype))
    if write_filtered:
        filtered_fh.close()

    return subs, indels

def insertSNPs(subs, indels, ref):
    if subs:
        new_chrs = insertSubs(subs, ref)
    if indels:
        new_chrs = insertIndels(indels, new_chrs)    
    return new_chrs

def insertSubs(subs, ref):
    ''' subs are 1 based, ref is 0'''
    new_chrs = dict()
    for chr in ref:
        new_chrs[chr] = list(ref[chr])

        chr_snp = [s for s in subs if s[0] == chr]
        #print len(chr_snp), chr
        for sub in chr_snp:
            pos = sub[1]
            alt = sub[3]
            genotype = sub[4]
            if '1/1' not in genotype:
                continue
            new_chrs[chr][pos - 1] = alt
        #print 'pre subs:', ref[chr]
        #print 'postsubs:', ''.join(new_chrs[chr])
    return new_chrs

def insertIndels(indels, ref):
    
    # construct a new chr adding insertions and ommiting deleted parts
    # if a deletion, make a new part from the old idx to pos, setting the new
    # pos to after the deletion
    # if an insertion, add a new part from idx to pos + the insertion and
    # setting the idx to after the insertion
    new_chrs = dict()
    for chr, seq in ref.items():
        chr_indels = [i for i in indels if i[0] == chr]
        chr_indels.sort()
        if not chr_indels:
            new_chrs[chr] = ''.join(seq)
            #print chr+' has no indels'
            continue

        p = [] 
        idx = 0
        net = 0
        for indel in chr_indels:
            pos = indel[1]
            ref = indel[2]
            alt = indel[3]
            
            if ',' in alt:
                continue #heterozygous
            #print 'before processing ',indel, p, idx
            if len(alt) < len(ref):
                #deletion
                del_len = len(ref)-1
                p.extend(seq[idx:pos])
                idx = pos + del_len
                net -= del_len
            elif len(alt) > len(ref):
                #insertion
                p.extend(seq[idx:pos+1]) # insertions goes after pos
                p.extend(list(alt[1:]))
                idx = pos
                net += len(alt)-1
            else:
                print 'error: alt has same len as ref',indel
                sys.exit()
            #print 'after:', p, idx
        # add last bit
        p.extend(seq[idx:])
        #print 'end:', p
        new_chrs[chr] = ''.join(p)
        #print len(chr_indels), chr, 'net change:', net
    #print 'post indels:', new_chrs
    return new_chrs

def writeFasta(chrs,out_file):
    fh = open(out_file, 'w')
    keys = chrs.keys()
    keys.sort()
    for chr in keys:
        fh.write('>'+chr+'\n')
        fh.write(chrs[chr]+'\n')
    fh.close()
    
def main(argv):
    if len(argv) < 2:
        print 'usage: script.py ref_seq.fasta sample -|vcf_file [filtered_vcf_out]'
        sys.exit()
    ref_file = argv[1]
    sample = argv[2]
    vcf_file = argv[3]
    if vcf_file == '-':
        fh = sys.stdin
    else:
        fh = open(vcf_file)
    out_vcf = None
    try:
        out_vcf = sys.argv[4]
    except:
        pass

    subs, indels = parseVCF(fh,write_filtered=out_vcf)
    print 'subs', len(subs)
    print 'indels', len(indels)
    if len(indels) != 0:
        print len(subs), len(indels), float(len(subs))/len(indels)

    chrs = readFasta(ref_file)
    new_chrs = insertSNPs(subs, indels, chrs)

    for chr in new_chrs:
        print chr, 'diff:', len(new_chrs[chr])-len(chrs[chr])

    keys = new_chrs.keys()
    for k in keys:
        new_header = k + ' ' + sample
        new_chrs[new_header] = new_chrs[k][:]
        del new_chrs[k]
    print new_chrs.keys()
    writeFasta(new_chrs, sample+'.fasta')


if __name__ == '__main__':
    main(sys.argv)
