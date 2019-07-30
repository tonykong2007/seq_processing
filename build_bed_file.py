# -*- coding: utf-8 -*-
"""
This script takes a list of transcript names or rs numbers and generates a sorted BED file.

Usage:
    build_bed_file.py [TRANSCRIPT-NAME-1 or RS_1] [TRANSCRIPT-NAME-2 OR RS_2]
"""
import argparse
import re
import sys

from autoprimer import Gene, SNP, split_transcript_name

def generate_bed(identifiers, flanking=10):
    unsorted_output = list()
    custom_region_count = 0
    for identifier in identifiers:
        if re.match(r'rs[0-9]+', identifier):
            # use the SNP class instead
            g = SNP(identifier, version='GRCh37')
            if int(g.start) == int(g.end):
                unsorted_output.append(['chr' + g.chromosome,
                                       g.start, g.start,
                                       g.name
                                       ]
                                       )
            else:
                for pos in range(int(g.start), int(g.end) + 1):
                    unsorted_output.append(['chr' + g.chromosome,
                                           pos, pos,
                                           g.name
                                           ]
                                           )
            continue
        elif re.match(r'chr[0-9]+:[0-9]+', identifier):
            custom_region_count += 1
            chromosome, pos = identifier.split(':')
            unsorted_output.append([chromosome,
                                    pos, pos,
                                    'custom_region_' + str(custom_region_count).zfill(3)
                                    ])
            continue
        elif re.search(r'-[0-9][0-9][0-9]', identifier):
            # user-specified transcript
            gene_name, transcript_number = split_transcript_name(identifier)
            g = Gene(gene_name, version='GRCh37')
            g.set_transcript(identifier)
        else:
            # use canonical transcript from Ensembl
            gene_name = identifier
            g = Gene(gene_name, version='GRCh37')
            g.set_transcript()
        exons = g.list_exon_regions()
        for i, exon in enumerate(exons, 1):
            cds = g.exon_to_translated(exon)
            if cds[1] is not None and cds[2] is not None:
                assert cds[1] < cds[2], 'cds[1] is not less than cds[2]!'
                flanked_cds = [cds[1] - flanking, cds[2] + flanking]
                unsorted_output.append(['chr' + str(cds[0]),
                                        flanked_cds[0], flanked_cds[1],
                                        gene_name.upper() + '.EX{exon:03}'.format(exon=i) +
                                        '.' + str(cds[1]) + '.' + str(cds[2])
                                        ]
                                        )
    return unsorted_output

def get_key(region):
    """Return the unique key generated from chromosome and start position"""
    return region[0].replace('chr', '').zfill(2) + '-' + str(region[1]).zfill(10)

def sort_regions(region_list):
    """Return the sorted region list"""
    print('# Fixing and sorting intervals...', file=sys.stderr)
    # fixes the potential issue of end > start in some indels
    for region in region_list: 
        if int(region[2]) < int(region[1]):
            region[1], region[2] = region[2], region[1]
    # actual sorting
    sorted_regions = sorted(region_list, key=lambda x: get_key(x))
    return sorted_regions

def dedup_regions(region_list):
    """Return the deduplicated region list"""
    print('# Deduplicating intervals...', file=sys.stderr)
    dedup_regions = dict()
    for region in region_list:
        if get_key(region) in dedup_regions:
            dedup_regions[get_key(region)][3] += ';' + region[3]
        else:
            dedup_regions[get_key(region)] = region
    dedup_region_list = list()
    for key in dedup_regions:
        dedup_region_list.append(dedup_regions[key])
    return dedup_region_list

def main():
    parser = argparse.ArgumentParser(description='Builds a BED file from a list of provided genes/ transcripts.')
    parser.add_argument('transcripts', nargs='+', help='gene name or trancript, e.g. ATP7B or ATP7B-201')
    parser.add_argument('-f', dest='cds_flank', type=int, default=10, help='flanking region to add to each target interval (default: 10)')
    args = parser.parse_args()
    # pass the command line arguments to the generate_bed() function
    unsorted_output = generate_bed(args.transcripts, flanking=args.cds_flank)
    dedup_output = dedup_regions(unsorted_output)
    sorted_output = sort_regions(dedup_output)
    for line in sorted_output:
        print('\t'.join([str(x) for x in line]))


if __name__ == '__main__':
    main()
