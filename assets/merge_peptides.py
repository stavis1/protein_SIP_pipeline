#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 15:51:33 2025

@author: 4vt
"""

import os
import dill
import pandas as pd

dillfiles = [f for f in os.listdir() if f.endswith('.dill')]
peptides = []
for dillfile in dillfiles:
    with open(dillfile, 'rb') as dillhandle:
        peptides.extend(dill.load(dillhandle))
with open('peptides.dill', 'wb') as dillhandle:
    dill.dump(peptides, dillhandle)

def read_tables(file):
    try:
        pd.read_csv(file, sep = '\t')
    except pd.errors.EmptyDataError:
        return None
tables = [read_tables(f) for f in os.listdir() if f.endswith('.tsv')]
tables = [t for t in tables if t is not None]
if tables:
    pd.concat(tables).to_csv('peptides.tsv', sep = '\t', index = False)
else:
    print('No valid enriched peptides found.')
    exit(1)
