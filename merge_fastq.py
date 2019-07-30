#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat 1 Jun 2019

Script for merging multiple gzipped FASTQ files

@author: HCC604
"""
import re
import sys
import os
import glob

# process the directory and diectory name
source_folder = input('Please drag-and-drop the folder containing the .fq.gz files, then press Enter:\n')
source_folder = source_folder.rstrip()
source_folder = source_folder.lstrip("'").rstrip("'")

print('# Now checking', source_folder, file=sys.stderr)
assert os.path.isdir(source_folder), 'Wrong input detected. Must be a folder.'

folder_name = os.path.basename(source_folder)
print('# Output files will be named:', folder_name + '[1|2]_.fq.gz', file=sys.stderr)

output_path = os.path.realpath(source_folder)
print('# Output would be written to:', output_path, file=sys.stderr)

# now prepare the file lists
forward_reads = glob.glob(os.path.join(output_path, '*_1.fq.gz'))
forward_reads.sort()
reverse_reads = glob.glob(os.path.join(output_path, '*_2.fq.gz'))
reverse_reads.sort()

# now merge the forward reads
print('# Now merging the forward reads...', file=sys.stderr)
forward_read_file_string = ''
for f in forward_reads:
    print('\t' + f + '...', file=sys.stderr)
    forward_read_file_string += '"' + f + '" '
merge_command_1 = 'gunzip -c ' + forward_read_file_string + '| gzip -c > ' + os.path.join(output_path, folder_name + '_1.fq.gz')
print(merge_command_1, file=sys.stderr)
os.system(merge_command_1)

# now merge the reverse reads
print('# Now merging the reverse reads...', file=sys.stderr)
reverse_read_file_string = ''
for r in reverse_reads:
    print('\t' + r + '...', file=sys.stderr)
    reverse_read_file_string += '"' + r + '" '
merge_command_2 = 'gunzip -c ' + reverse_read_file_string + '| gzip -c > ' + os.path.join(output_path, folder_name + '_2.fq.gz')
print(merge_command_2, file=sys.stderr)
os.system(merge_command_2)

print('# Done!', file=sys.stderr)


