#! /usr/bin/env python
import os
import sys
import sequenceTools as st

def parseSNPs(snp_file):
    subs = list()
    indels = list()
    #strain = snp_file.split('.')[0]

    for line in snp_file.readlines():
        if line.startwith('#'):
            if line.startswith('#CHROM'):
                strain = line.split()[9].split('.')[0]
            continue
        p = line.split('\t')
        chr =p[0]
        pos = int(p[1])-1 # pileup is 1-based
        ref = p[3]
        var = p[4]
        if ref == '*':
            indels.append((chr,pos, var))
        else:
            subs.append((chr, pos, var))
    return subs, indels

def insertSNPs(subs, indels, ref):
    new_chrs = insertSubs(subs, ref)
    new_chrs = insertINDELs(indels, new_chrs)    
    return new_chrs

def insertSubs(subs, ref):
    new_chrs = dict()
    for chr in ref:
        new_chrs[chr] = list(ref[chr])

        chr_snp = [s for s in subs if s[0] == chr]
        print len(chr_snp), chr
        for s in chr_snp:
            new_chrs[chr][s[1]] = s[2]
        #print 'pre subs:', ref[chr]
        #print 'postsubs:', ''.join(new_chrs[chr])
    return new_chrs

def insertINDELs(indels, ref):
    
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
            print chr+' has no indels'
            continue

        p = [] 
        idx = 0
        net = 0
        for indel in chr_indels:
            a,b = indel[2].split('/')
            if a != b:
                continue #heterozygous
            #print 'before processing ',indel, p, idx
            if indel[2].startswith('-'):
                #deletion
                del_len = len(a)-1
                p.extend(seq[idx:indel[1]])
                idx = indel[1] + del_len
                net -= del_len
            else:
                #insertion
                p.extend(seq[idx:indel[1]+1]) # insertions goes after pos
                p.extend(list(a[1:]))
                idx = indel[1]+1
                net += len(a)-1
            #print 'after:', p, idx
        # add last bit
        p.extend(seq[idx:])
        #print 'end:', p
        new_chrs[chr] = ''.join(p)
        print len(chr_indels), chr, 'net change:', net
    #print 'post indels:', new_chrs
    return new_chrs

def writeFasta(chrs,out_file):
    fh = open(out_file, 'w')
    for chr in chrs:
        fh.write('>'+chr+'\n')
        fh.write(chrs[chr]+'\n')
    fh.close()
    
def main(argv):
    if len(argv) < 2:
        print 'usage: script.py ref_seq.fasta -|vcf_file'
        sys.exit()
    snp_file = argv[1]
    if len(argv) == 4:
        fh = open(vcf_file)
    else:
        fh = sys.stdin
    
    subs, indels, strain = parseSNPs(fh)
    print 'subs', len(subs)
    print 'indels', len(indels)
    if len(indels) != 0:
        print len(subs), len(indels), float(len(subs))/len(indels)

    chrs = st.readFasta(ref_file)
    print chrs.keys()
    for k in chrs:
        header = k.split()
        arm = header[0]
        new_header = strain+'_'+arm
        chrs[new_header] = chrs[k]
        del chrs[k]
    print chrs.keys()
    
    new_chrs = insertSNPs(subs, indels, chrs)
    for chr in new_chrs:
        print chr, 'diff:', len(new_chrs[chr])-len(chrs[chr])
    writeFasta(new_chrs, 'new.fasta')


if __name__ == '__main__':
    main(sys.argv)
