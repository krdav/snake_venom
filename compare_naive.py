#!/usr/bin/env python
from __future__ import division, print_function
import sys, os, glob, csv, random, copy, time, shutil, re
import matplotlib
# matplotlib.use('Agg')
matplotlib.use('PDF')
# random.seed(666)
from itertools import cycle
csv.field_size_limit(sys.maxsize)
from anarci import anarci
from collections import Counter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
import numpy as np
partis_path = '/fh/fast/matsen_e/kdavidse/partis_master/partis'
sys.path.insert(1, partis_path + '/python')
import utils
import glutils
sys.path.insert(1, '/fh/fast/matsen_e/kdavidse/aammp')
import pandas as pd
import seaborn as sns
from matplotlib import pyplot
import jellyfish
try:
    def hamming_distance(s1, s2):
        if s1 == s2:
            return 0
        else:
            return jellyfish.hamming_distance(s1, s2)
    assert(hamming_distance('ABC', 'ABCD') == 1)
except:
    def hamming_distance(s1, s2):
        if s1 == s2:
            return 0
        else:
            return jellyfish.hamming_distance(unicode(s1), unicode(s2))
    assert(hamming_distance('ABC', 'ABCD') == 1)

### Notice a gap is added as the 21th amino acid:
AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
AA_INDEX = {aa:i for i, aa in enumerate(AA_LIST)}

AHO_L = 149


class BadNaive(Exception):
    '''When the naive sequence is bad.'''


class CustomCry(Exception):
    '''When the naive sequence is bad.'''


def repair_seq(seq, naiveDNA, vj_bounds, keep_check=False):
    '''
    Repair unobserved sequence with the naive sequence, codon by codon,
    then strip the Ns if the sequence is still N padded and return the in-frame sequence.
    This function also checks the validity of a sequence if the keep_check option is set true.
    '''

    naiveDNA = naiveDNA[vj_bounds[0]:vj_bounds[1]]
    if len(naiveDNA)%3 != 0:
        naiveDNA = naiveDNA[0:-(len(naiveDNA)%3)]  # Remove the untranslated end

    seq = seq[vj_bounds[0]:vj_bounds[1]]
    if len(seq)%3 != 0:
        seq = seq[0:-(len(seq)%3)]  # Remove the untranslated end
    assert(len(seq) == len(naiveDNA))

    # Convert to mutable:
    naiveDNA = list(naiveDNA)
    trim_seq = list(seq)
    # Repair codon by codon, starting by the 3-prime end.
    # Reverse the codons to the the 3-prime end:
    naiveDNA = naiveDNA[::-1]
    trim_seq = trim_seq[::-1]
    # Repair until no more consecutive Ns from N padding:
    for i in range(3, MAX_REP, 3):
        if 'N' in trim_seq[i-3:i]:
            trim_seq[i-3:i] = naiveDNA[i-3:i]
        else:
            break
    # Back reverse again, then trim the 5-prime end the same way:
    naiveDNA = naiveDNA[::-1]
    trim_seq = trim_seq[::-1]
    for i in range(3, MAX_REP, 3):
        if 'N' in trim_seq[i-3:i]:
            trim_seq[i-3:i] = naiveDNA[i-3:i]
        else:
            break
    # Join to string and check for inconsistencies:
    trim_seq = ''.join(trim_seq)
    # If in "check for keeping" mode, return either true or false:
    if keep_check is True:
        if 'N' in trim_seq or '*' in str(Seq(trim_seq, generic_dna).translate()) or 3*MIN_LEN > len(trim_seq):
            # print('Fifth skip')
            return False
        else:
            return True
    # Else expect no inconsistencies and return trimmed sequence if none:
    # by using [:-9] stop codon is allowed in the last 3 codons
    if 'N' in trim_seq or '*' in str(Seq(trim_seq[:-9], generic_dna).translate()) or 3*MIN_LEN > len(trim_seq):
        raise BadNaive('Problem with the naive sequence. Either contains N or stop codon, or is too short: {}'.format(trim_seq))
    else:
        return trim_seq


def repair_seq_debug(seq, naiveDNA, vj_bounds, keep_check=True):
    '''
    Repair unobserved sequence with the naive sequence, codon by codon,
    then strip the Ns if the sequence is still N padded and return the in-frame sequence.
    This function also checks the validity of a sequence if the keep_check option is set true.
    '''

    naiveDNA = naiveDNA[vj_bounds[0]:vj_bounds[1]]
    if len(naiveDNA)%3 != 0:
        naiveDNA = naiveDNA[0:-(len(naiveDNA)%3)]  # Remove the untranslated end

    seq = seq[vj_bounds[0]:vj_bounds[1]]
    if len(seq)%3 != 0:
        seq = seq[0:-(len(seq)%3)]  # Remove the untranslated end
    assert(len(seq) == len(naiveDNA))

    # Convert to mutable:
    naiveDNA = list(naiveDNA)
    trim_seq = list(seq)
    # Repair codon by codon, starting by the 3-prime end.
    # Reverse the codons to the the 3-prime end:
    naiveDNA = naiveDNA[::-1]
    trim_seq = trim_seq[::-1]
    # Repair until no more consecutive Ns from N padding:
    for i in range(3, MAX_REP, 3):
        if 'N' in trim_seq[i-3:i]:
            trim_seq[i-3:i] = naiveDNA[i-3:i]
        else:
            break
    # Back reverse again, then trim the 5-prime end the same way:
    naiveDNA = naiveDNA[::-1]
    trim_seq = trim_seq[::-1]
    for i in range(3, MAX_REP, 3):
        if 'N' in trim_seq[i-3:i]:
            trim_seq[i-3:i] = naiveDNA[i-3:i]
        else:
            break
    # Join to string and check for inconsistencies:
    trim_seq = ''.join(trim_seq)
    # If in "check for keeping" mode, return either true or false:
    if keep_check is True:
        if 'N' in trim_seq or '*' in str(Seq(trim_seq, generic_dna).translate()) or 3*MIN_LEN > len(trim_seq):
            # print('Fifth skip')
            print(trim_seq)
            return False
        else:
            return True
    # Else expect no inconsistencies and return trimmed sequence if none:
    # by using [:-9] stop codon is allowed in the last 3 codons
    if 'N' in trim_seq or '*' in str(Seq(trim_seq[:-9], generic_dna).translate()) or 3*MIN_LEN > len(trim_seq):
        raise BadNaive('Problem with the naive sequence. Either contains N or stop codon, or is too short: {}'.format(trim_seq))
    else:
        return trim_seq


def Most_Common(lst):
    '''Return the most used list element.'''
    data = Counter(lst)
    return data.most_common(1)[0][0]


def extract_seqs(fnam, uid2iso):
    '''Reads a partis cluster-annotations files and extrats relevant information and sequences.'''
    # Read cluster annotations into a data list of dictionaries:
    with open(fnam) as fh:
        reader = csv.DictReader(fh)
        data = list(reader)

    sequences_i = list()
    info_i = list()
    for row in data:
        fnam_base = fnam.split('_partitions')[0]
        cwd = os.getcwd()
        if 'IgK' in fnam_base:
            locus = 'igk'
        elif 'IgL' in fnam_base:
            locus = 'igl'
        else:
            locus = 'igh'
        # Process the partis data row and add germline information:
        try:
            utils.process_input_line(row)
            # Read default germline info
            glfo = glutils.read_glfo('{}/_output/{}/hmm/germline-sets'.format(cwd, fnam_base), locus=locus)
            utils.add_implicit_info(glfo, row)
        except Exception as e:  # Skip rows that cannot be processed
            print('First skip')
            print(e)
            continue

        uids = [dl + [u] if (len(dl) > 0 and dl[0] != '') else [u] for dl, u in zip(row['duplicates'], row['unique_ids'])]

        # Extract the full N padded naive sequence,
        # and find the v -and j gene bound on this naive sequence:
        cdr3_bounds = (row['codon_positions']['v'], row['codon_positions']['j'] + 3)
        vj_bounds = (row['regional_bounds']['v'][0], row['regional_bounds']['j'][1])
        if row['invalid'] is True or (cdr3_bounds[0]-cdr3_bounds[1])%3 != 0:
            print('Invalid clonal family, skipping.')
            continue

        naiveDNA = row['naive_seq']
        if repair_seq(naiveDNA, naiveDNA, vj_bounds, keep_check=True) is False:  # Skip naive sequences too short or with stop codons:
            # print('Third skip')
            if len(row['input_seqs'][:]) > 100:
                print('Bad naive even after 100 seqs in clonal family.')
                repair_seq_debug(naiveDNA, naiveDNA, vj_bounds)
            continue
        trimmed_naiveDNA = repair_seq(naiveDNA[:], naiveDNA[:], vj_bounds)
        naiveAA = str(Seq(trimmed_naiveDNA, generic_dna).translate())

        # There has been a name change and this try/except
        # is meant to provide backwards compatability:
        try:
            lseq = row['input_seqs'][:]
        except:
            lseq = row['seqs'][:]
        ir_lseq = row['indel_reversed_seqs']
        stop_seq = row['stops']
        assert(len(lseq) == len(ir_lseq))
        assert(len(lseq) == len(stop_seq))
        # Only keep sequences without indels and stop codons and minimum length amino acid length:
        ### ir_lseq[i] == '' or lseq[i] == ir_lseq[i]  <-- No indels
        ### stop_seq[i]  <-- No partis annotated stops (there seems still to be stops after these are removed though)
        ### repair_seq(lseq[i], naiveDNA, vj_bounds, keep_check=True)  <-- Checks whether the sequence is long enougth or have stop codons
        keep_idx = [1 if ((ir_lseq[i] == '' or lseq[i] == ir_lseq[i]) and stop_seq[i] is False and repair_seq(lseq[i], naiveDNA, vj_bounds, keep_check=True)) else 0 for i in range(len(lseq))]

        # Now only keep those sequences that passed QC:
        lseq = [s for s, keep in zip(lseq, keep_idx) if keep == 1]
        # Exclude small clonal families:
        if len(lseq) < MIN_OBS:
            # print(len(lseq))
            # print('Fourth skip')
            continue
        # Get amino acid sequences:
        lAAseq = [str(Seq(repair_seq(s[:], naiveDNA[:], vj_bounds), generic_dna).translate()) for s in lseq]
#        mut_freqs = [s for s, keep in zip(row['mut_freqs'], keep_idx) if keep == 1]
#        print(row['n_mutations'].split(':'))
        Nmuts = [int(s) for s, keep in zip(row['n_mutations'].split(':'), keep_idx) if keep == 1]
        abundance = [len(d) for d, keep in zip(uids, keep_idx) if keep == 1]
        uids = [s for s, keep in zip(uids, keep_idx) if keep == 1]
        assert(len(Nmuts) == len(lseq))
        assert(len(abundance) == len(lseq))
        assert(len(uids) == len(lseq))
#        assert(len(mut_freqs) == len(lseq))
        # Convert frequency to counts and throw out info for discarded sequences:
#        Nmuts = [int(round(float(t[0])*len(t[1].strip('N')))) for i, t in enumerate(zip(mut_freqs, lseq))]

        # Deduplicate AAseqs and lseq according to the AA deduplication:
        '''
        lAAseq_dict = dict()
        lAAseq_sort = dict()
        lseq_dedup = list()
        for i, aa in enumerate(lAAseq):
            if aa in lAAseq_sort:
                lAAseq_sort[aa].append((i, repair_seq(lseq[i][:], naiveDNA[:], vj_bounds), abundance[i]))
            else:
                lAAseq_sort[aa] = [(i, repair_seq(lseq[i][:], naiveDNA[:], vj_bounds), abundance[i])]

        for i, aa in enumerate(lAAseq_sort):
            lAAseq_dict[aa] = [t[0] for t in lAAseq_sort[aa]]
            s = sorted(lAAseq_sort[aa], )
            ab_seq = sorted(lAAseq_sort[aa], key=lambda x: x[2], reverse=True)[0][1]
            lseq_dedup.append(ab_seq)

        assert(len(lAAseq_dict) == len(lseq_dedup))
        # Make the deduplicated list and take the mutation rates,
        #  as the mutation rate for the deduplicated sequence:
        lAAseq_dedup = list()
        Nmuts_dedup = list()
        abundance_dedup = list()
        for aa, idxs in lAAseq_dict.items():
            lAAseq_dedup.append(aa)
            Nmut_list = [float(Nmuts[i]) for i in idxs]
            Nmuts_dedup.append(int(round(sum(Nmut_list)/len(Nmut_list))))
            abundance_list = [abundance[i] for i in idxs]
            abundance_dedup.append(sum(abundance_list))
        assert(len(lAAseq_dedup) == len(Nmuts_dedup))
        assert(len(lAAseq_dedup) == len(abundance_dedup))
        assert(len(lAAseq_dedup) == len(lseq_dedup))

        # Exclude small clonal families:
        if len(lAAseq_dedup) < MIN_OBS:
            # print(len(lseq))
            # print('Fourth skip')
            continue
        '''
        iso_list = [[uid2iso[u] for u in ul] for ul in uids]
        # Store the results in a list:
        sequences_i.append(['naive_seq', naiveAA])  # This format is for ANARCI numbering
        info_i.append({'fnam': fnam, 'v_gene': row['v_gene'], 'd_gene': row['d_gene'], 'j_gene': row['j_gene'],
                       'naive_seq': naiveAA, 'naive_seq_DNA': trimmed_naiveDNA, 'Nmuts': Nmuts[:], 'abundance': abundance[:],
                       'AAseqs': lAAseq[:], 'DNAseqs': lseq[:], 'UID': uids[:], 'isotype': iso_list[:],
                       'CDR3_start': cdr3_bounds[0], 'CDR3_end': cdr3_bounds[1]})

    return(sequences_i, info_i)


def allele2grp(s):
    if 'IGLDx' in s or 'IGKDx' in s:
        return s
    try:
        m = re.match('IG[KLH][VDJ][0-9]+', s)
        return m.group(0)
    except:
        print(s)
        sys.exit()


def make_dataframe(info_i_j):
    header = ['clusterID', 'naiveAA', 'vdj_len', 'naive', 'Nseqs', 'CDR3_start', 'CDR3_end', 'v_grp', 'd_grp', 'j_grp', 'v_gene', 'd_gene', 'j_gene', 'sample', 'locus', 'filename', 'input_seqs', 'input_seqsAA', 'Nmuts', 'abundance', 'UID', 'isotype']
    df = [header]
    for i, info in enumerate(info_i_j):
        if 'IgK' in info['fnam']:
            locus = 'IgK'
        elif 'IgL' in info['fnam']:
            locus = 'IgL'
        else:
            locus = 'IgH'

        cols = [i, info['naive_seq'], len(info['naive_seq']), info['naive_seq_DNA'], len(info['Nmuts']), info['CDR3_start'], info['CDR3_end'],  allele2grp(info['v_gene']), allele2grp(info['d_gene']), allele2grp(info['j_gene']), info['v_gene'], info['d_gene'], info['j_gene'], info['fnam'].split('_')[0], locus, info['fnam'], ':'.join(info['DNAseqs']), ':'.join(info['AAseqs']), ':'.join(list(map(str,info['Nmuts']))), ':'.join(list(map(str, info['abundance']))), ':'.join(['@'.join(uidl) for uidl in info['UID']]), ':'.join(['@'.join(iso) for iso in info['isotype']])]
        try:
            assert(len(cols) == len(header))
        except:
            print(cols)
            print(header)
            print(len(cols))
            print(len(header))
            raise Exception('More data columns than header elements.')
        df.append(cols)
    return df


def write_dataframe(df, outfile):
    fh_out = open(outfile, 'w')
    for row in df:
        fh_out.write(','.join(list(map(str, row))))
        fh_out.write('\n')
    fh_out.close()


def run_file(package):
    '''
    Wrapping the sequence extraction etc. for each cluster file.
    This function can then be called outside by a pool of subprocesses.
    '''
    f, uid2iso = package
    sequences_i, info_i = extract_seqs(f, uid2iso)
    assert(len(sequences_i) == len(info_i))
    print('Reading {} containing {} GCs.'.format(f, len(sequences_i)))
    if len(sequences_i) == 0:
        return False
    return info_i


def AUC(res, locus, st, cut):
    counts = Counter(res[locus][st])
    hd = [k for k, v in counts.items() if k <= cut]
    c = [v for k, v in counts.items() if k <= cut]
    hd, c = zip(*sorted(zip(hd, c)))
    acc = [sum(c[0:i]) + c[i] for i in range(len(c))]
    sc = sum(v for k, v in counts.items())
    auc = sum(acc[i]*(hd[i+1]-hd[i]) for i in range(len(hd)-1)) / sc
    return(auc)


def naiveAA_plot(res, xlimit, trim, downsample):
#    %matplotlib inline
    from matplotlib.backends.backend_pdf import PdfPages
    ds = '_downsampled' if downsample else ''
    pp = PdfPages('NaiveAA_comparison_trim{}{}.pdf'.format(trim, ds))
    for locus in res:
        fig, ax = pyplot.subplots(figsize=(14,10))
        for st in res[locus]:
            if st[0][0:3] == '3FT' and st[1][0:3] == '3FT':
                color = 'green'
            elif st[0][0:3] == 'PLA' and st[1][0:3] == 'PLA':
                color = 'blue'
            else:
                color = 'red'
            auc = AUC(res, locus, st, xlimit)
            ax = sns.distplot(
                res[locus][st],
                label='{} vs. {} AUC {:.1f}'.format(st[0], st[1], auc),
                color=color,
                kde=False,
                bins=list(range(0, 1000)),
                norm_hist=True,  # On/off to normalize y-axis
                hist_kws={'histtype':'step', 'cumulative':True, 'lw':3}
            )
            lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            ax.set_title('Locus: {}'.format(locus))
            ax.set_xlim(0, xlimit)
        fig.savefig(pp, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    pp.close()


def compare_naive(table, trim=0, outfile_p=False, downsample=False):
    if outfile_p is not False:
        headers = table.pop(0)
        df = pd.DataFrame(table, columns=headers)
        df = df.loc[:, ['Nseqs', 'locus', 'sample', 'v_grp', 'd_grp', 'j_grp', 'vdj_len', 'naiveAA']]
        df.to_pickle('{}_pickle'.format(outfile_p))
    else:
        df = table.copy(deep=True)

    # Trim out small clonal families:
    df = df[df['Nseqs'] > trim]

    row_part = dict()
    dist_table = dict()

    samples = list(set(df['sample']))
    loci = set(df['locus'])

    if downsample:
        for locus in loci:
            smallest = 9999999999
            for sample in samples:
                print(len(df[(df['locus'] == locus) & (df['sample'] == sample)]))
                size = sum(np.array(df['locus'] == locus) & np.array(df['sample'] == sample))
                if size < smallest:
                    smallest = size
            print('Downsampling to:', smallest)
            for sample in samples:
                # Finding the downsample:
                ds = df[(df['locus'] == locus) & (df['sample'] == sample)].sample(smallest, replace=False)
                # Dropping all columns:
                df = df[np.invert(np.array(((df['locus'] == locus) & (df['sample'] == sample))))]
                # Adding the downsample:
                df = df.append(ds)
    v_grps = set(df['v_grp'])
    d_grps = set(df['d_grp'])
    j_grps = set(df['j_grp'])
    vdj_len = set(df['vdj_len'])

            
    for index, row in df.iterrows():
        # kt = (row['locus'], row['v_grp'], row['d_grp'], row['j_grp'], row['vdj_len'])
        kt = (row['locus'], row['v_grp'], row['d_grp'], row['j_grp'])
        if kt not in row_part:
            row_part[kt] = [index]
        else:
            row_part[kt].append(index)
    sumstat = [len(l) for l in row_part.values()]
#    print('These are the number of entries in each partition:', list(map(str, sumstat)))

    # res[locus][s1:s2] = [dist, 5, 2, 0,...]
    res = dict()
    for locus in loci:
        res[locus] = dict()
        for i in range(len(samples)):
            for j in range(i+1, len(samples)):
                si = samples[i]
                sj = samples[j]
                res[locus][(si, sj)] = list()
                for kt, li in row_part.items():
                    if locus != kt[0]:
                        continue
                    naivei = list(df.loc[li][df.loc[li]['sample'] == si]['naiveAA'])
                    naivej = list(df.loc[li][df.loc[li]['sample'] == sj]['naiveAA'])
                    if len(naivei) == 0 or len(naivej) == 0:
                        continue
                    hdi = [min(hamming_distance(ni, nj) for ni in naivei) for nj in naivej]
                    hdj = [min(hamming_distance(ni, nj) for nj in naivej) for ni in naivei]
                    res[locus][(si, sj)].extend(hdi)
                    res[locus][(si, sj)].extend(hdj)
    naiveAA_plot(res, 30, trim, downsample)


def get_uid2isotype(fnams):
    uid2isotype = dict()
    for fnam in fnams:
        with open(fnam) as fh:
            d = dict()
            for l in fh:
                 if l.startswith('>'):
                     uid = l[1:].split('|')[0]
                     m = re.search('-Ig\w+-', l)
                     iso = m.group(0).strip('-')
                     if uid not in d:
                         d[uid] = iso
                     else:
                         print('UID multiple occurence in file:', uid)
        sample = fnam.split('_')[0]
        uid2isotype[sample] = d
    return uid2isotype


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Create a BCR sequences profile from partis clonal family partitions.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--glob_pat', type=str, required=True, default='*cluster-annotations.csv', help='Glob pattern to find the partis "cluster annotation" files.')
    parser.add_argument('--glob_pat_iso', type=str, required=True, default=None, help='Isotype info.')
    parser.add_argument('--nsubs', type=int, required=False, default=200, help='How many subsampled profiles should be taken out?')
    parser.add_argument('--fsub', type=float, required=False, default=None, help='Fraction of sequences in the full profile to include in the subsample.')
    parser.add_argument('--totN_sub', type=int, required=False, default=10, help='Total number of sequences in the full profile to include in the subsample.')
    parser.add_argument('--nproc', type=int, required=False, default=1, help='Number of processes to start.')
    parser.add_argument('--dataset_name', type=str, required=False, default='not_specified', help='Dataset name.')
    parser.add_argument('--outfile_p', type=str, required=False, default='DATASET_NAMEtopXXX_aammp_profiles.txt', help='Output name for profile.')
    parser.add_argument('--outfile_s', type=str, required=False, default='DATASET_NAMEtopXXX_aammp_profiles_seqs.txt', help='Output name for seqeunces belonging to each profile.')
    parser.add_argument('--outfile_sub_p', type=str, required=False, default='DATASET_NAMEtopXXX_aammp_profiles_subYYY.txt', help='Output name for subsampled profile.')
    parser.add_argument('--outfile_sub_s', type=str, required=False, default='DATASET_NAMEtopXXX_aammp_profiles_seqs_subYYY.txt', help='Output name for sequences belonging to each subsampled profile.')
    parser.add_argument('--MIN_OBS', type=int, required=False, default=1, help='Minimum number of sequences in each clonal family.')
    parser.add_argument('--MIN_LEN', type=int, required=False, default=70, help='Minimum length of any given amino acid sequence in the dataset.')
    parser.add_argument('--MAX_REP', type=int, required=False, default=30, help='Max number of nt. repaired at each end.')
    global args
    args = parser.parse_args()

    # Make cmd arg later:
    global SIM_SIZE
    SIM_SIZE = 10

    global MIN_OBS
    global MIN_LEN
    global MAX_REP
    MIN_OBS = args.MIN_OBS
    MIN_LEN = args.MIN_LEN
    MAX_REP = args.MAX_REP
    nsubs = args.nsubs
    fsub = args.fsub
    totN_sub = args.totN_sub
    assert((fsub, totN_sub).count(None) == 1)

    global SPECIES
    SPECIES = 'mouse'

    # Insert dataset name in filename:
    outfile_p = args.outfile_p.replace('DATASET_NAME', '{}'.format(args.dataset_name))
    outfile_s = args.outfile_s.replace('DATASET_NAME', '{}'.format(args.dataset_name))
    outfile_sub_p = args.outfile_sub_p.replace('DATASET_NAME', '{}'.format(args.dataset_name))
    outfile_sub_s = args.outfile_sub_s.replace('DATASET_NAME', '{}'.format(args.dataset_name))
    # Insert testset size in filename:
    outfile_p = outfile_p.replace('topXXX', '{}'.format(str(nsubs)))
    outfile_s = outfile_s.replace('topXXX', '{}'.format(str(nsubs)))
    outfile_sub_p = outfile_sub_p.replace('topXXX', '{}'.format(str(nsubs)))
    outfile_sub_s = outfile_sub_s.replace('topXXX', '{}'.format(str(nsubs)))
    # Insert subsample size or fraction in filename:
    if totN_sub is None:
        outfile_p = outfile_p.replace('YYY', '{}'.format(str(fsub)))
        outfile_s = outfile_s.replace('YYY', '{}'.format(str(fsub)))
        outfile_sub_p = outfile_sub_p.replace('YYY', '{}'.format(str(fsub)))
        outfile_sub_s = outfile_sub_s.replace('YYY', '{}'.format(str(fsub)))
    else:
        outfile_p = outfile_p.replace('YYY', 'N{}'.format(str(totN_sub)))
        outfile_s = outfile_s.replace('YYY', 'N{}'.format(str(totN_sub)))
        outfile_sub_p = outfile_sub_p.replace('YYY', 'N{}'.format(str(totN_sub)))
        outfile_sub_s = outfile_sub_s.replace('YYY', 'N{}'.format(str(totN_sub)))

    glob_res = list(glob.glob(args.glob_pat))
    print('These are the files that are going to be parsed:\n{}'.format('\n'.join(glob_res)))

    glob_isores = list(glob.glob(args.glob_pat_iso))
    print('These are the files that are going to be parsed for isotype info:\n{}'.format('\n'.join(glob_isores)))
    uid2isotype = get_uid2isotype(glob_isores)

    packages = [(f, uid2isotype[f.split('_')[0]]) for f in glob_res]
    if args.nproc == 1:
        results = map(run_file, packages)  # Run without subprocess
    else:
        import multiprocessing
        # Paralellize the process:
        pool = multiprocessing.Pool(args.nproc)  # Start the pool
        results = pool.map_async(run_file, packages, chunksize=1)  # Run subprocesses
        pool.close()
        pool.join()
        results = results.get()

    # Unpack and merge the results:
    info_i_j = list()           # GC data for j in file i extra information beloning to the sequences
    for t in results:
        if t is not False:
            info_i_j.extend(t)

    print('Total GCs in files:', len(info_i_j))
    print('Total sequences in all GCs:', sum([len(ll['AAseqs']) for ll in info_i_j]))

    # Sort lists according to number of sequence in each clonal family:
    sort_idx = [t[0] for t in sorted(enumerate(info_i_j), key=lambda x: len(x[1]['AAseqs']), reverse=True)]
    info_i_j = [info_i_j[i] for i in sort_idx]

    # Generate dataframes with the profiles and extra info (notice the slice from nsubs excluded):
    df = make_dataframe(info_i_j)

    # Write the profiles:
    write_dataframe(df, outfile_p)

    compare_naive(copy.deepcopy(df), trim=0, outfile_p=outfile_p, downsample=True)
    compare_naive(copy.deepcopy(df), trim=5, outfile_p=outfile_p, downsample=True)
    compare_naive(copy.deepcopy(df), trim=0, outfile_p=outfile_p, downsample=False)
    compare_naive(copy.deepcopy(df), trim=5, outfile_p=outfile_p, downsample=False)


if __name__ == '__main__':
    pretty_random_fnam = str(random.randint(1, 10**100))
    global TMPDIR
    TMPDIR = '/tmp/kd_tmp_' + pretty_random_fnam
    os.mkdir(TMPDIR)  # Make a tmp dir to dump crap
    main()
    shutil.rmtree(TMPDIR)  # rm -rf tmp dir
    print('Done')

#    try:
#        main()
#        shutil.rmtree(TMPDIR)  # rm -rf tmp dir
#    except Exception as e:
#        shutil.rmtree(TMPDIR)  # rm -rf tmp dir
#        print('There was an error:', e)


# ANARCI --sequence some_fasta.fa --outfile some_fasta.anno --scheme aho --restrict H --ncpu 1 --use_species human
# unique_ids,v_gene,d_gene,j_gene,cdr3_length,mut_freqs,input_seqs,indel_reversed_seqs,naive_seq,indelfos,duplicates,v_per_gene_support,d_per_gene_support,j_per_gene_support,v_3p_del,d_5p_del,d_3p_del,j_5p_del,v_5p_del,j_3p_del,vd_insertion,dj_insertion,fv_insertion,jf_insertion,mutated_invariants,in_frames,stops
# ((117, ' '), 'V')
# unique_ids,v_gene,d_gene,j_gene,cdr3_length,mut_freqs,input_seqs,indel_reversed_seqs,naive_seq,indelfos,duplicates,v_per_gene_support,d_per_gene_support,j_per_gene_support,v_3p_del,d_5p_del,d_3p_del,j_5p_del,v_5p_del,j_3p_del,vd_insertion,dj_insertion,fv_insertion,jf_insertion,mutated_invariants,in_frames,stops

# partis
# /fh/fast/matsen_e/kdavidse/partis/bin/partis run-viterbi --chain h --species human --infname inp.fasta --outfname out.csv
# cmd = '{}/bin/partis partition --chain {} --species {} --infname tmp/{}.fa --outfname tmp/{}.csv'.format(partis_path, chain, species, inpf, outf)
# out-cluster-annotations.csv
# out.csv
