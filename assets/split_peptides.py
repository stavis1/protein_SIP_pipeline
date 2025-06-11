#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 15:45:15 2025

@author: 4vt
"""

import sys
import dill

dillfile = sys.argv[1]
Nchunks = int(sys.argv[2])
with open(dillfile, 'rb') as dillhandle:
    step, peptides = dill.load(dillhandle)

subsets = [[] for _ in range(Nchunks)]
for i, pep in enumerate(peptides):
    subsets[i%Nchunks].append(pep)

for i, subset in enumerate(subsets):
    with open(f'subset_step{step}_{i}.dill', 'wb') as dillhandle:
        dill.dump((step, subset), dillhandle)