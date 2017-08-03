#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
IgA
IgD
IgE
IgGA
IgGb
IgK
IgL
IgM
'''


import glob
import re, random
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import MultipleSeqAlignment

gp = '*.fasta'
split_size = 1000000
isotypes = ['IgK', 'IgL', 'IgH']
# iso_dict = {iso:[] for iso in isotypes}


def dump_split(basename, sp, iso, seqlist):
    fho = open('{}_{}_split_{}.fa'.format(basename, iso, sp), 'w')
    print 'Writing:', '{}_{}_split_{}.fa'.format(basename, iso, sp)
    for i, s_obj in enumerate(seqlist):
        # fho.write('>{}\n{}\n'.format(i, str(s_obj.seq)))
        fho.write('>{}\n{}\n'.format(s_obj.id, str(s_obj.seq)))
    fho.close()


for f in glob.glob(gp):
    print 'Processing:', f
    basename = f[:-6]

    seqs = list(SeqIO.parse(f, 'fasta'))
    # random.shuffle(seqs)  # In place shuffle

    i = 1
    sp = 1
    Nseqs = len(seqs)
    iso_dict = {iso:[] for iso in isotypes}
    iso = ''
    for s in seqs:
        dnaseq = str(s.seq)
        dnaseq = dnaseq.upper()
        dnaseq = dnaseq.strip('N')
        if re.match('^[ATGCN]+$', dnaseq):  #### Added N to keep all sequences
            if 'IgK' in  s.id:
                iso = 'IgK'
            elif 'IgL' in  s.id:
                iso = 'IgL'
            else:
                iso = 'IgH'
            # Notice here how seqlist is used like a pointer back to the dict list:
            seqlist = iso_dict[iso]
            seqlist.append(s)
            i += 1
            if len(seqlist) % split_size == 0:
                dump_split(basename, sp, iso, seqlist)
                seqlist[:] = list()
                sp += 1
        else:
            pass
            #print dnaseq
    for iso in isotypes:
       seqlist = iso_dict[iso]
       dump_split(basename, sp, iso, seqlist)
    print 'Starting with:', Nseqs
    print 'Ending with:', i-1



