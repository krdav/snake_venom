#! /usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import re
import sys, os


gp = '*.fa'

qsub_temp = 'sub_temp.qsub'
#node_spec = '32:fatnode'
node_spec = '32:fatnode'
nproc = '32'

with open(qsub_temp) as fh:
    qsub_string = fh.read()

locus = ''
for f in glob.glob(gp):
    basename = f[:-3]
    qsub_out = basename+'.qsub'
    qsub_string_cp = qsub_string[:]
    if 'IgK' in  basename:
        locus = 'igk'
    elif 'IgL' in  basename:
        locus = 'igl'
    else:
        locus = 'igh'
    with open(qsub_out, 'w') as fho:
        qsub_string_cp = qsub_string_cp.replace('@BASENAME@', basename)
        qsub_string_cp = qsub_string_cp.replace('@fatORthinNODE@', node_spec)
        qsub_string_cp = qsub_string_cp.replace('@NPROC@', nproc)
        qsub_string_cp = qsub_string_cp.replace('@LOCUS@', locus)
        fho.write(qsub_string_cp)
    os.system('qsub {}'.format(qsub_out))

