#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

psm_file = next(f for f in os.listdir() if f.endswith('.psm.txt'))
data = pd.read_csv(psm_file, sep = '\t', comment = '#')
data['file'] = [f.replace('.FT2', '.mzML') for f in data['Filename']]
data['seq'] = [s[1:-1] for s in data['IdentifiedPeptide']]
data['proteins'] = [p[1:-1] for p in data['ProteinNames']]
data = data[data['TargetMatch'] == 'T']
data.to_csv(sys.argv[1] + '.psms', sep = '\t', index = False)
