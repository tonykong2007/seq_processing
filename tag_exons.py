#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 15:52:50 2017

@author: HCC604
"""
import math
import os
import re
import sys

import numpy as np

from Bio import SeqIO, pairwise2
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

from autoprimer import Gene, get_flanking_regions, get_sequence
from web_search import save_mutation_taster_pdf

CDS_FLANK = 10
SCORE_THRESHOLD = 50
SUPPRESS_VARIANT_THRESHOLD = 5

def parse_ABI(path):
    record = SeqIO.read(path, 'abi')
    G = record.annotations['abif_raw']['DATA9']
    A = record.annotations['abif_raw']['DATA10']
    T = record.annotations['abif_raw']['DATA11']
    C = record.annotations['abif_raw']['DATA12']
    RUNT2 = record.annotations['abif_raw']['RUNT2']
    RUND2 = record.annotations['abif_raw']['RUND2']
    PLOC2 = record.annotations['abif_raw']['PLOC2']
    PBAS2 = record.annotations['abif_raw']['PBAS2']
    P2BA1 = record.annotations['abif_raw']['P2BA1']
    SMPL1 = record.annotations['abif_raw']['SMPL1']
    PCON1 = record.annotations['abif_raw']['PCON1']

    record_dict = {
            'A': A,
            'T': T,
            'C': C,
            'G': G,
            'peaks': PLOC2,
            'bases': [x for x in PBAS2],
            'quality': [ord(x) for x in PCON1],
            'bases2': P2BA1,
            'sample': SMPL1,
            'date': RUND2,
            'time': RUNT2
            }
    return record_dict

def reverse_complement(input_abi_dict):
    if 'bases_anno' in input_abi_dict:
        raise ValueError('Reverse complement of annotated abi_dict is not supported.')
    rc_abi_dict = dict()
    rc_abi_dict['A'], rc_abi_dict['T'] = input_abi_dict['T'][::-1], input_abi_dict['A'][::-1]
    rc_abi_dict['G'], rc_abi_dict['C'] = input_abi_dict['C'][::-1], input_abi_dict['G'][::-1]
    last_pos = len(input_abi_dict['A'])
    rc_abi_dict['peaks'] = [last_pos - x for x in input_abi_dict['peaks']]
    rc_abi_dict['peaks'] = rc_abi_dict['peaks'][::-1]
    trans = str.maketrans('ATGCYRSWKMBDHVN', 'TACGRYSWMKVHDBN')
    rc_abi_dict['bases'] = [x.translate(trans) for x in input_abi_dict['bases']]
    rc_abi_dict['bases'] = rc_abi_dict['bases'][::-1]
    rc_abi_dict['bases2'] = [x.translate(trans) for x in input_abi_dict['bases2']]
    rc_abi_dict['bases2'] = rc_abi_dict['bases2'][::-1]
    rc_abi_dict['quality'] = input_abi_dict['quality'][::-1]
    return rc_abi_dict

def align_ABI(abi_dict, cds_seqs, exon_cds_start):
    ''' perform local alignment to find matching segments '''
    abi_seq = ''.join(abi_dict['bases'])
    per_exon_alignment = list()
    for i, cds_seq in enumerate(cds_seqs, 1):
        if cds_seq is None:
            per_exon_alignment.append({'exon': i,
                                   'abi': '',
                                   'cds': '',
                                   'score': 0,
                                   'cds_start': exon_cds_start[i-1]
                                   })
            continue
        alignments = pairwise2.align.localms(abi_seq, cds_seq, 1, -2, -5, -2)
        abi_aligned = alignments[0][0]
        cds_aligned = alignments[0][1]
        alignment_score = alignments[0][2]
        per_exon_alignment.append({'exon': i,
                                   'abi': abi_aligned,
                                   'cds': cds_aligned,
                                   'score': alignment_score,
                                   'cds_start': exon_cds_start[i-1]
                                   })
    return per_exon_alignment

def generate_annotation_symbols(actual_sequence, reference_sequence, cds_flank,
                                exon_number, exon_cds_start, cds_tick_interval=5):
    symbols = list()
    non_gap_length = len(str(reference_sequence).replace('-', ''))
    head_flank_start = 1
    head_flank_end = cds_flank
    print('# -' + str(cds_flank) + ':', head_flank_start, head_flank_end, file=sys.stderr)
    tail_flank_start = non_gap_length - cds_flank + 1
    tail_flank_end = non_gap_length
    print('# -' + str(cds_flank) + ':', tail_flank_start, tail_flank_end, file=sys.stderr)
    this_cds_base_count = 0
    cds_count = exon_cds_start
    if len(actual_sequence) != len(reference_sequence):
        raise ValueError('Aligned sequences not equal in length!')
    for actual_base, reference_base in zip(actual_sequence, reference_sequence):
        if actual_base == '-' and reference_base == '-':
            # possible alignment error
            raise ValueError('Alignment error. Check pairwise alignment output.')
        if actual_base == '-' and reference_base != '-':
            # possible deletion
            this_cds_base_count += 1
            # this is for handling a CDS region deletion (as it is not handled later)
            if head_flank_end < this_cds_base_count < tail_flank_start:
                cds_count += 1
            symbols.append('d')
            print('::DEBUG::')
            print(this_cds_base_count, cds_count, actual_base, reference_base)
        if reference_base == '-' and actual_base != '-':
            # either a base not covered by reference, or an insertion
            symbols.append('.')
        if actual_base != '-' and reference_base != '-':
            # a successfully aligned base
            this_cds_base_count += 1
            if head_flank_start <= this_cds_base_count <= head_flank_end:
                if actual_base != str(reference_base).upper():
                    symbols.append('!')
                else:
                    symbols.append('-')
            elif tail_flank_start <= this_cds_base_count <= tail_flank_end:
                if actual_base != str(reference_base).upper():
                    symbols.append('!')
                else:
                    symbols.append('+')
            else:
                # inside the CDS region of the exon
                cds_count += 1 # this is 1-based; minus 1 when printing
                if actual_base != str(reference_base).upper():
                    if (cds_count - 1) == 1 or (cds_count - 1) % cds_tick_interval == 0:
                        symbols.append('!' + ' c.' + str(cds_count - 1))
                    else:
                        symbols.append('!')
                else:
                    if (cds_count - 1) == 1 or (cds_count - 1) % cds_tick_interval == 0:
                        symbols.append('E' + str(exon_number) + ' c.' + str(cds_count - 1))
                    else:
                        symbols.append('E' + str(exon_number))

    print('>aligned|actual_sequence')
    print(actual_sequence)
    print('>aligned|reference_sequence')
    print(reference_sequence)
    print(len(actual_sequence), len(symbols), file=sys.stderr)
    return symbols

def resolve_ambiguous_base(sample_base, reference_base):
    '''Resolve ambiguous bases into mutated bases'''
    resolved_base = None
    print(reference_base, '>', sample_base, file=sys.stderr)
    sample_base = str.upper(sample_base)
    reference_base = str.upper(reference_base)
    if (sample_base == 'A' or sample_base == 'T'
        or sample_base == 'C' or sample_base == 'G'):
        resolved_base = sample_base
    elif sample_base == 'R':
        if reference_base == 'A':
            resolved_base = 'G'
        elif reference_base == 'G':
            resolved_base = 'A'
        else:
            print('# sample_base coerced to A')
            resolved_base = 'A'
    elif sample_base == 'Y':
        if reference_base == 'C':
            resolved_base = 'T'
        elif reference_base == 'T':
            resolved_base = 'C'
        else:
            print('# sample_base coerced to C')
            resolved_base = 'C'
    elif sample_base == 'S':
        if reference_base == 'C':
            resolved_base = 'G'
        elif reference_base == 'G':
            resolved_base = 'C'
        else:
            print('# sample_base coerced to G')
            resolved_base = 'G'
    elif sample_base == 'W':
        if reference_base == 'A':
            resolved_base = 'T'
        elif reference_base == 'T':
            resolved_base = 'A'
        else:
            print('# sample_base coerced to T')
            resolved_base = 'T'
    elif sample_base == 'M': # bug fix on 20190327, previously missing implementation
        if reference_base == 'C':
            resolved_base = 'A'
        elif reference_base == 'A':
            resolved_base = 'C'
        else:
            print('# sample_base coerced to C')
            resolved_base = 'C'
    elif sample_base == 'K': # bug fix on 20190327, previously missing implementation
        if reference_base == 'T':
            resolved_base = 'G'
        elif reference_base == 'G':
            resolved_base = 'T'
        else:
            print('# sample_base coerced to T')
            resolved_base = 'T'
    else:
        print('Warning: mutatant base', sample_base, ' will be coerced to non-reference base', file=sys.stderr)
        if reference_base == 'A':
            resolved_base = 'G'
        elif reference_base == 'T':
            resolved_base = 'C'
        elif reference_base == 'C':
            resolved_base = 'T'
        elif reference_base == 'G':
            resolved_base = 'A'
        else:
            raise ValueError('Reference base', reference_base, 'invalid!')
    return resolved_base

def find_variants(aligned_sample, aligned_reference):
    if len(aligned_sample) != len(aligned_reference):
        raise ValueError('Alignment error!')
    # find indels
    indels = list()
    dels = re.finditer('[A-Z]+([-]+)[A-Z]+', str.upper(aligned_sample))
    for d in dels:
        del_start = d.start(1)
        del_end = d.end(1)
        if del_start >= 20:
            indels.append(aligned_reference[del_start-20:del_start] + '[' +
                          aligned_reference[del_start:del_end] + '/' +
                          '-' + ']' + aligned_reference[del_end:del_end+20])
        else:
            indels.append('[' +
                          aligned_reference[del_start:del_end] + '/' +
                          '-' + ']' + aligned_reference[del_end:del_end+20])
    ins = re.finditer('[A-Z]+([-]+)[A-Z]+', str.upper(aligned_reference))
    for i in ins:
        ins_start = i.start(1)
        ins_end = i.end(1)
        if ins_start >= 20:
            indels.append(aligned_reference[ins_start-20:ins_start] + '[' +
                          '-' + '/' +
                          aligned_sample[ins_start:ins_end] +
                          ']' + aligned_reference[ins_end:ins_end+20])
        else:
            indels.append('[' +
                          '-' + '/' +
                          aligned_sample[ins_start:ins_end] +
                          ']' + aligned_reference[ins_end:ins_end+20])
    # find SNPs
    snps = list()
    for i, j in enumerate(aligned_sample):
        k = aligned_reference[i]
        if j != '-' and k != '-':
            if str.upper(j) != str.upper(k):
                if i >= 20:
                    snps.append(aligned_reference[i-20:i] + '[' + k + '/' + resolve_ambiguous_base(j, k) + ']' + aligned_reference[i+1:i+21])
                else:
                    snps.append('[' + k + '/' + resolve_ambiguous_base(j, k) + ']' + aligned_reference[i+1:i+21])
    variants = snps + indels
    if len(variants) > 0:
        for v in variants:
            print('MT <= Variant snippet:', v, file=sys.stderr)
    return variants

def annotate_ABI(abi_dict, per_exon_alignment, score_threshold=SCORE_THRESHOLD,
                 tick_interval=5):
    anno_abi_dict = dict(abi_dict)
    variants = list()
    anno_abi_dict['bases_anno'] = [str(i+1) for i in range(len(abi_dict['bases']))]
    for i, exon in enumerate(per_exon_alignment):
        if exon['score'] >= score_threshold:
            print('# Alignment score of exon', i+1, 'is', exon['score'], file=sys.stderr)
            print('# Annotation will be performed', file=sys.stderr)
            print("# abi_dict['bases']=", len(abi_dict['bases']), file=sys.stderr)
            print("# exon['abi'] non gap =", len(str(exon['abi']).replace('-', '')), file=sys.stderr)
            annotation_symbols = generate_annotation_symbols(exon['abi'], exon['cds'], CDS_FLANK,
                                                             exon['exon'], exon['cds_start'])
            base_count = 0
            for j, base in enumerate(exon['abi']):
                if base == '-':
                    continue
                else:
                    base_count += 1
                    try:
                        if int(anno_abi_dict['bases_anno'][base_count - 1]) % tick_interval != 0:
                            anno_abi_dict['bases_anno'][base_count - 1] = ''
                    except:
                        pass
                    anno_abi_dict['bases_anno'][base_count - 1] += ' ' + annotation_symbols[j]
            variants = list(variants) + find_variants(exon['abi'], exon['cds'])
    return (anno_abi_dict, variants)

def cut_ABI(input_abi_dict, per_row=60):
    '''Cut abi_dict into rows, each with per_row bases'''
    cutoffs = list()
    for i, peak in enumerate(input_abi_dict['peaks']):
        if i % per_row == 0:
            cutoffs.append({'base': i, 'peak_loc': peak})
    rows = list()
    peak_lrange = 0
    base_lrange = 0
    for i in range(1, len(cutoffs)):
        this_row = dict()
        peak_urange = cutoffs[i]['peak_loc']
        base_urange = cutoffs[i]['base']
        for base in ['A', 'T', 'C', 'G']:
            this_row[base] = input_abi_dict[base][peak_lrange:peak_urange]
        this_row['peaks'] = input_abi_dict['peaks'][base_lrange:base_urange]
        this_row['bases'] = input_abi_dict['bases'][base_lrange:base_urange]
        if 'bases_anno' in input_abi_dict:
            this_row['bases_anno'] = input_abi_dict['bases_anno'][base_lrange:base_urange]
        rows.append(this_row)
        peak_lrange = peak_urange
        base_lrange = base_urange
    # final row
    this_row = dict()
    for base in ['A', 'T', 'C', 'G']:
        this_row[base] = input_abi_dict[base][peak_urange:]
    this_row['peaks'] = input_abi_dict['peaks'][base_urange:]
    this_row['bases'] = input_abi_dict['bases'][base_urange:]
    if 'bases_anno' in input_abi_dict:
            this_row['bases_anno'] = input_abi_dict['bases_anno'][base_urange:]
    rows.append(this_row)
    return rows

def percentile_height(abi_dict, percentile=99):
    '''Return peak height at peak_pos'''
    A_height = abi_dict['A']
    T_height = abi_dict['T']
    C_height = abi_dict['C']
    G_height = abi_dict['G']
    height = A_height + T_height + C_height + G_height
    return np.percentile(height, percentile)

def plot_ABI(plt_object, abi_dict_chunk, chunk_index):
    plt_object.plot(abi_dict_chunk['A'], color='green', linewidth=0.5)
    plt_object.plot(abi_dict_chunk['T'], color='red', linewidth=0.5)
    plt_object.plot(abi_dict_chunk['C'], color='blue', linewidth=0.5)
    plt_object.plot(abi_dict_chunk['G'], color='black', linewidth=0.5)
    if chunk_index == 0:
        offset = 0
    else:
        offset = abi_dict_chunk['peaks'][0]
    # annotate the peaks
    for i, peak_pos in enumerate([x - offset for x in abi_dict_chunk['peaks']]):
        base = abi_dict_chunk['bases'][i]
        base_anno = abi_dict_chunk['bases_anno'][i]
        if base == 'A':
            color = 'green'
        elif base == 'T':
            color =  'red'
        elif base == 'C':
            color = 'blue'
        elif base == 'G':
            color = 'black'
        else:
            color = 'pink'
        # background colouring based on annotation
        border_colour = 'none'
        base_bg = 'white'
        if 'E' in base_anno:
            base_bg = 'yellow'
        if '-' in base_anno or '+' in base_anno:
            base_bg = 'cyan'
        if '!' in base_anno:
            border_colour = 'red'
        plt_object.text(peak_pos,
                        percentile_height(abi_dict_chunk),
                        base,
                        horizontalalignment='center',
                        color=color,
                        name='Courier New',
                        fontweight='bold',
                        bbox={'facecolor':base_bg, 'alpha':0.5, 'pad':0.1, 'edgecolor':border_colour})
    plt.sca(plt_object)
    plt.xticks([x - offset for x in abi_dict_chunk['peaks']],
               abi_dict_chunk['bases_anno'],
               rotation='vertical',
               fontsize=8,
               name='Calibri')
    return (abi_dict_chunk['peaks'], abi_dict_chunk['bases'], [x - offset for x in abi_dict_chunk['peaks']], abi_dict_chunk['bases_anno'])

def multi_plot_ABI(abi_dict_chunks, pdf_object,
                   abi_dict, is_reverse_complement,
                   gene_name, transcript,
                   source_file, rows_per_page=7):
    plotted = list()
    not_saved = None
    if is_reverse_complement:
        strand = 'REVERSE STRAND'
    else:
        strand = 'FORWARD STRAND'
    total_pages = math.ceil(len(abi_dict_chunks) / rows_per_page)
    this_page = 0
    print('Total', len(abi_dict_chunks), 'chunks', file=sys.stderr)
    for i, chunk in enumerate(abi_dict_chunks):
        print('now plotting chunk', i, file=sys.stderr)
        # open a new page
        if len(plotted) % rows_per_page == 0:
            if (len(abi_dict_chunks) - len(plotted)) >= rows_per_page:
                rows_this_page = rows_per_page
            else:
                rows_this_page = (len(abi_dict_chunks) - len(plotted)) % rows_per_page
            if rows_this_page == 1:
                rows_this_page = 2 # the minimum number for AxesSubplot to work
            print('...plotted thus far:', len(plotted), file=sys.stderr)
            print('Creating a new page with', rows_this_page, 'rows', file=sys.stderr)
            f, axarr = plt.subplots(rows_this_page, 1, figsize=(8, rows_this_page * 11/7), dpi=300)
            this_page += 1
            not_saved = True
        plotted.append(plot_ABI(axarr[i % rows_per_page], chunk, i))
        plt.tight_layout(rect=[0, 0.03, 1, 0.9])
        for item in ['sample', 'date', 'time']:
            if item not in abi_dict:
                abi_dict[item] = 'undefiled'
                print('Warning: ', item, 'not found in abi_dict!', file=sys.stderr)
        plt.suptitle('Sample: ' + abi_dict['sample'] + ' '
                     'Date: ' + abi_dict['date'] + ' ' +
                     'Time: ' + abi_dict['time'] + '\n' +
                     'Filename: ' + os.path.basename(sys.argv[1]) + '\n' +
                     strand + ' ' +
                     'Gene: ' + gene_name + ' ' +
                     'Transcript: ' + transcript + ' '
                     'Page: ' + str(this_page) + ' of ' + str(total_pages),
                     fontweight='bold',
                     name='Arial',
                     fontsize=12)
        # closes the page if rows_per_page rows have been filled
        if len(plotted) > 0 and len(plotted) % rows_per_page == 0:
            pdf_object.savefig()
            not_saved = False
    # save the figure even if the page was not completely filled
    if not_saved:
        pdf_object.savefig()
    return plotted

if __name__ == '__main__':
    abi_path = sys.argv[1]
    # obtain coding +- 10 bp sequences for chromatogram annotation
    print(sys.argv, file=sys.stderr)
    gene_name, transcript_no = sys.argv[2].split('-')    
    this_gene = Gene(gene_name)
    print(this_gene.list_transcripts(), file=sys.stderr)
    assert transcript_no.isdigit() == True
    transcript = gene_name + '-' + transcript_no
    print(this_gene.set_transcript(transcript), file=sys.stderr)
    exons = this_gene.list_exon_regions()
    coding_regions = [this_gene.exon_to_translated(er) for er in exons]
    cds_seqs = list()
    cds_pos_count = 1
    exon_cds_start = list() # this is the lowest c. number in the exon
    for cds in coding_regions:
        if cds[1] is None or cds[2] is None:
            cds_seqs.append(None)
            exon_cds_start.append(cds_pos_count)
            continue
        upstream, downstream = get_flanking_regions(cds, flank=CDS_FLANK)
        coreseq = cds
        flanked_cds_seq = (get_sequence(upstream) +
                           get_sequence(coreseq) +
                           get_sequence(downstream))
        cds_seqs.append(flanked_cds_seq)
        exon_cds_start.append(cds_pos_count)
        cds_pos_count += len(get_sequence(coreseq))

    print('# Now converting file...', file=sys.stderr)
    abi_dict = parse_ABI(abi_path)
    print('# Now performing alignment for forward strand...', file=sys.stderr)
    alignment_dict = align_ABI(abi_dict, cds_seqs, exon_cds_start)
    print('# Reverse complementation...', file=sys.stderr)
    abi_dict_rc = reverse_complement(abi_dict)
    print('# Now performing alignment for reverse strand...', file=sys.stderr)
    alignment_dict_rc = align_ABI(abi_dict_rc, cds_seqs, exon_cds_start)

    annotated_ABI = None
    rc_status = None
    # automatically determine of directionality
    print('# Current annotation alignment score threshold:', SCORE_THRESHOLD, file=sys.stderr)
    forward_score = sum([x['score'] for x in alignment_dict])
    print('# Forward strand score:', forward_score, file=sys.stderr)
    reverse_score = sum([x['score'] for x in alignment_dict_rc])
    print('# Reverse strand score:', reverse_score, file=sys.stderr)
    if forward_score > SCORE_THRESHOLD and forward_score > reverse_score:
        print('# FORWARD strand selected.', file=sys.stderr)
        rc_status = False
        annotated_ABI, variants = annotate_ABI(abi_dict, alignment_dict)
    elif reverse_score > SCORE_THRESHOLD and reverse_score > forward_score:
        print('# REVERSE strand selected.', file=sys.stderr)
        rc_status = True
        annotated_ABI, variants = annotate_ABI(abi_dict_rc, alignment_dict_rc)
    else:
        print('# Cannot determine which strand to use. Quitting...', file=sys.stderr)
        sys.exit()

    print('# Annotation done. Now plotting...', file=sys.stderr)

    rows = cut_ABI(annotated_ABI)

    with PdfPages(os.path.basename(sys.argv[1]) + '.pdf') as pdf:
        p = multi_plot_ABI(rows, pdf, annotated_ABI, rc_status,
                           gene_name, transcript, sys.argv[1])
    print('Chromatogram PDF saved!', file=sys.stderr)

    try:
        # MutationTaster gene name, if provided
        gene_name = sys.argv[4]
    except:
        print('Using gene_name:', gene_name, file=sys.stderr)

    if len(variants) <= SUPPRESS_VARIANT_THRESHOLD:
        for i, variant in enumerate(variants, 1):
            save_mutation_taster_pdf(gene_name, this_gene.transcript_id,
                                     variant,
                                     os.path.basename(sys.argv[1]) + '.variant-' + format(i, '03'),
                                     os.path.basename(sys.argv[1]) + '.variant-' + format(i, '03') + '.pdf')
    else:
        print('# Too many variants, skipping...', file=sys.stderr)



