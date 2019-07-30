#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 11:18:19 2017

@author: Tom
"""

import sys

import os.path

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from scipy import misc

def parse_vcf(vcf_path):
    """Return a list of PNG filenames and corresponding VCF lines, with a caller tag"""
    png_list = list()
    is_freebayes = False
    is_bcftools = False
    with open(vcf_path, 'r') as f:
        variant_count = 0
        for line in f:
            if str(line).startswith('#'):
                if str(line).startswith('##source=freeBayes'):
                    is_freebayes = True
                elif str(line).startswith('##bcftoolsVersion'):
                    is_bcftools = True
            else:
                variant_count += 1
                line = str(line).split('\t')
                chrom = line[0]
                pos = int(line[1])
                varid = line[2]
                ref = line[3]
                alt = line[4]
                qual = float(line[5])
                filtertag = line[6]
                info = line[7]
                info_dict = dict()
                for item in info.split(';'):
                    if item == 'ALLELE_END':
                        pass
                    else:
                        key, value = item.split('=')
                    info_dict[key] = value
                if 'SAMPLE' in info_dict:
                    # read sample name from INFO field
                    sample_name = info_dict['SAMPLE']
                else:
                    # infer sample name from VCF filename
                    sample_name = os.path.basename(vcf_path)
                    sample_name = str(sample_name).split('.')[0]
                    info_dict['SAMPLE'] = sample_name
                png_filename = sample_name + 's' + str(variant_count) + '__' + chrom + '_' + str(pos) + '.png'
                png_list.append({
                                'filename': png_filename,
                                'chrom': chrom,
                                'pos': pos,
                                'varid': varid,
                                'ref': ref,
                                'alt': alt,
                                'qual': qual,
                                'filter': filtertag,
                                'info': info_dict
                                })
    return (png_list, is_freebayes, is_bcftools)

def short(string, l=15):
    """"""
    if len(string) > l:
        return string[0:l-1]
    else:
        return string

def w(info_dict, key):
    """Special function for printing INFO field values"""
    if key in info_dict:
        return info_dict[key]
    else:
        return 'N.A.'

def dp4(dp4_field, count):
    """Special function to return ref/ alt count from bcftools output"""
    ref_f, ref_r, alt_f, alt_r = dp4_field.split(',')
    if count == 'ref':
        return int(ref_f) + int(ref_r)
    elif count == 'alt':
        return int(alt_f) + int(alt_r)
    else:
        raise NotImplementedError

def parse_annotation(forcecall_regions):
    """Simple function to convert the annotation file to a list"""
    anno = list()
    with open(forcecall_regions, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                anno.append(line.rstrip().split('\t'))
    return anno


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print('Usage: ngs_variant_report.py [in.vcf] [forcecall.regions]')
        sys.exit()

    vcf_path = os.path.realpath(sys.argv[1])
    vcf_folder = os.path.dirname(vcf_path)
    png_list, called_by_fb, called_by_bcftools = parse_vcf(vcf_path)
    if len(sys.argv) == 3:
        rs_anno = parse_annotation(sys.argv[2])
        assert len(rs_anno) == len(png_list), "Annotation file mismatch!"

    # create a multi-page PDF file object for writing
    with PdfPages(os.path.join(vcf_path + '.pdf')) as pdf:
        variant_count = 0
        for png in png_list:
            variant_count += 1
            print('Reporting variant', variant_count, file=sys.stderr)
            png_filename = png['filename']
            png_path = os.path.join(vcf_folder + '/snapshots', png_filename) # the PNG files must be stored in a "snapshots" sub-directory
            this_image = misc.imread(png_path)
            fig = plt.figure(figsize=(11.69, 8.27), dpi=150)
            a = fig.add_subplot(1, 2, 1)
            a.set_title('IGV plot of ' + short(png['info']['SAMPLE']) + ' ' + png['chrom'] + ':' + str(png['pos']))
            plt.imshow(this_image)
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            a = fig.add_subplot(1, 2, 2)
            a.set_title('Variant summary')
            a.text(0.1, 0.95, 'Chromosome: ' + png['chrom'])
            a.text(0.1, 0.90, 'Position: ' + str(png['pos']))
            a.text(0.1, 0.85, 'Variant: ' + png['ref'] + '>' + png['alt'])
            if called_by_fb:
                a.text(0.1, 0.80, 'Frequency: ' + png['info']['AF']
                        + ' (' + png['info']['AO']
                        + '/' + png['info']['DP'] +
                        ')'
                        )
            elif called_by_bcftools:
                a.text(0.1, 0.80, 'Frequency: ' 
                        + ' ref=' + str(dp4(png['info']['DP4'], count='ref'))
                        + ' alt=' + str(dp4(png['info']['DP4'], count='alt')) 
                        + ' (total=' + str(dp4(png['info']['DP4'], count='ref') + dp4(png['info']['DP4'], count='alt')) + ')'
                        )
                a.text(0.1, 0.70, 'NCBI rs number (dbSNP150): ' + rs_anno[variant_count-1][3]) # the annotation list is 0-based
                a.text(0.1, 0.65, ' '.join(rs_anno[variant_count-1][0:3]))
            else:
                a.text(0.1, 0.80, 'Frequency: ' + png['info']['AF']
                        + ' (' + png['info']['VD']
                        + ' out of ' + png['info']['DP'] +
                        ' reads)'
                        )
                a.text(0.1, 0.75, 'Variant quality score: ' + str(png['qual']) + ' [' + png['filter'] + ']')
                a.text(0.1, 0.70, 'NCBI rs number (dbSNP150): ' + png['info']['avsnp150'])
                a.text(0.1, 0.65, 'ExAC: ' + png['info']['ExAC_ALL'] + ' [ALL] '
                                       + png['info']['ExAC_EAS'] + ' [EAS] ')
                a.text(0.1, 0.60, 'gnomAD exome: ' + png['info']['gnomAD_exome_ALL'] + ' [ALL] '
                                       + png['info']['gnomAD_exome_EAS'] + ' [EAS] ')
                a.text(0.1, 0.55, '1000G (2015Aug): ' + png['info']['1000g2015aug_all'] + ' [ALL] '
                                       + png['info']['1000g2015aug_eas'] + ' [EAS] ')
                # a.text(0.1, 0.50, 'Polyphen2 (HDIV): ' + png['info']['Polyphen2_HDIV_score']
                #                           + ' [' + png['info']['Polyphen2_HDIV_pred'] + ']')
                a.text(0.1, 0.45, 'SIFT: ' + png['info']['SIFT_score']
                                       + ' [' + png['info']['SIFT_pred'] + ']')
                a.text(0.1, 0.40, 'MutationTaster: ' + png['info']['MutationTaster_score']
                                       + ' [' + png['info']['MutationTaster_pred'] + ']')
                a.text(0.1, 0.35, 'InterVar: ' + png['info']['InterVar_automated'])
                a.text(0.1, 0.325, 'ClinVar: ' + png['info']['CLNSIG'] + '/'
                                           + png['info']['CLNDN'][0:30]) # truncate CLNDN at 30 characters
                # a.text(0.1, 0.30, png['info']['CLNACC'] + '/'
                #                    + png['info']['CLNDSDB'] + '/'
                #                    + png['info']['CLNDSDBID'])
                a.text(0.1, 0.275, 'refGene: ' + png['info']['Gene.refGene']
                                            + ' [' + png['info']['ExonicFunc.refGene'] + ']')
                a.text(0.1, 0.25, 'refGene AA change: ' + str(png['info']['AAChange.refGene']).replace(',', '\n'), va='top')
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()


