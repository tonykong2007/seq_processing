#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:56:35 2018

@author: HCC604
"""

import os
import sys

script_dir = os.path.dirname(os.path.realpath(__file__))

def set_proxy():
    proxy_txt_path = script_dir + '/proxy.txt'
    if os.path.isfile(proxy_txt_path):
        with open(proxy_txt_path, 'r') as f:
            for line in f:
                if line.startswith('set HTTP'):
                    os.system(line.rstrip())
                    print(line.rstrip(), file=sys.stderr)
