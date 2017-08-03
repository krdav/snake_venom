#! /usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import re
import sys, os


gp = '*-cluster-annotations.csv'
gawk_cmd = '''
gawk -F',' '{ if ($2!="") {print $0}}'
'''
gawk_cmd = gawk_cmd.strip()


for f in glob.glob(gp):
    basename = f[:-4]
    os.system('{} {} > {}_no_failed.csv'.format(gawk_cmd, f, basename))


