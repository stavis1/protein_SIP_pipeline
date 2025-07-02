#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 10:44:52 2025

@author: 4vt
"""

import sys
import re
import pandas as pd

config_file = sys.argv[1]

line_accum = []
flag = False
with open(config_file, 'r') as cfg:
    for line in cfg:
        if line.startswith('Element_List'):
            flag = True
        elif line.startswith('Element_Masses'):
            flag = False
        if flag:
            line_accum.append(line)

headers = re.findall(r'[A-Z][a-z]?(?=(?:,|\Z))', line_accum[0])
def residue(line):
    line = re.sub(r'\s*#.*\n?', '', line)
    if not line:
        return None
    res = re.search(r'\{([^\}])\}', line)
    if res is None:
        return None
    res = res.group(1)
    counts = re.findall(r'\d+(?=(?:,|\Z))', line)
    cols = {'code':res}
    cols.update({h:c for h,c in zip(headers, counts)})
    return pd.Series(cols)

residues = [residue(l) for l in line_accum[1:]]
pd.DataFrame([r for r in residues if r is not None]).to_csv('AA_formulae.tsv', sep = '\t', index = False)
