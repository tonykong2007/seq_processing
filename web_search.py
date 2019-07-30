#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:09:02 2017

@author: HCC604
"""
import pdfkit
import sys

from autoprimer import wget, check_url

def search_mutation_taster(gene, ensembl_transcript_id,
                           seq_snippet, alteration_name):
    '''Output HTML content of MutationTaster search result'''
    base_url = 'http://www.mutationtaster.org/cgi-bin/MutationTaster/MutationTaster69.cgi'
    query = base_url + '?'
    query += 'gene=' + gene + '&'
    query += 'transcript_stable_id_text=' + ensembl_transcript_id + '&'
    query += 'sequence_snippet=' + seq_snippet + '&'
    query += 'alteration_name=' + alteration_name
    # test the URL
    status_code = check_url(query)
    if status_code == 200:
        print('# Processing alteration ' + alteration_name, file=sys.stderr)
        return wget(query)
    else:
        return None

def save_mutation_taster_pdf(gene, ensembl_transcript_id,
                           seq_snippet, alteration_name,
                           filename):
    '''Perform MutationTaster search and save the webpage as PDF file'''
    result = search_mutation_taster(gene, ensembl_transcript_id,
                           seq_snippet, alteration_name)
    if result:
        pdfkit.from_string(result, filename)
        print('Mutation Taster result saved!', file=sys.stderr)
    else:
        print('# No valid result saved.', file=sys.stderr)