#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 14:19:37 2017

@author: HCC604
"""
import glob
import os
import sys

import setup_proxy
setup_proxy.set_proxy()

if len(sys.argv) != 3:
    print('usage: batch_tag_exons.py [PATH_TO_AB1_FOLDER] [GENE-TRANSCRIPT]', file=sys.stderr)
    print('e.g. batch_tag_exons.py c:/ATP7B/2018/patient1 ATP7B-201')
else:
    folder = sys.argv[1]
    ABI_files = glob.glob(folder+'/**/*.ab1', recursive=True)

    for file in ABI_files:
        print('# Now processing:', file, file=sys.stderr)
        os.system('python tag_exons.py "' + file + '" ' +
                  sys.argv[2])

