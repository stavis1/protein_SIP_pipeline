#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:36:34 2026

@author: stavis
"""

import os
import re
import subprocess
import numpy as np
import pandas as pd
from functools import cache
from collections import Counter

raw_data = pd.concat([pd.read_csv(f, sep = '\t', comment = '#') for f in os.listdir() if f.endswith('.sip')]).reset_index()
raw_data['SpecId'] = ['@'.join(str(f) for f in fields) for fields in zip(raw_data['Filename'],
                                                                         raw_data['SearchName'],
                                                                         raw_data['ScanNumber'],
                                                                         raw_data['ParentCharge'],
                                                                         raw_data['Rank'])]
data = raw_data.copy()

#metadata columns needed for percolator
data['Label'] = [-1 if 'Rev_' in p else 1 for p in data['ProteinNames']]
data['Peptide'] = data['IdentifiedPeptide']
data['Proteins'] = [('DECOY_' if 'Rev_' in p else '')+p for p in data['ProteinNames']]
data['ScanNr'] = data['ScanNumber']

#predictors for percolator
@cache
def getseq(p):
    return re.search(r'\[([^]]+)\]', p).group(1)
data['seq'] = [getseq(p) for p in data['IdentifiedPeptide']]
pred_cols = ['MeasuredParentMass',
             'scanNrCentroidDistance',
             'Rank',
             'Score',
             'peplen',
             'N_missed',
             'mass_err',
             'label_percent',
             'Nenz',
             'Cenz',
             'seqCount',
             'scanCount',
             'maxProtCount',
             'deltaScore',
             'chargeCount']
for charge in set(data['ParentCharge']):
    colname = f'charge_{charge}'
    pred_cols.append(colname)
    data[colname] = np.array(data['ParentCharge'] == charge, dtype = int)
data['peplen'] = [len(s) for s in data['seq']]
data['N_missed'] = [len(re.findall('[KR](?=.)', s)) for s in data['seq']]
data['mass_err'] = data['CalculatedParentMass'] - data['MeasuredParentMass']
data['label_percent'] = [int(re.search(r'_(\d+)Pct', n).group(1))/100000 for n in data['SearchName']]
data['Nenz'] = [int(s.startswith('[')) for s in data['IdentifiedPeptide']]
data['Cenz'] = [int(s.endswith(']')) for s in data['IdentifiedPeptide']]
seqCount = Counter(data['seq'])
scanCount = Counter(data['ScanNr'])
protCount = Counter(p for n in data['ProteinNames'] for p in n[1:-1].split(','))
data['seqCount'] = [np.log(seqCount[s]) for s in data['seq']]
data['scanCount'] = [np.log(scanCount[s]) for s in data['ScanNr']]
data['maxProtCount'] = [np.log(max(protCount[p] for p in n[1:-1].split(','))) for n in data['ProteinNames']]
maxscores = data.groupby('ScanNr')['Score'].apply(np.max).to_dict()
data['deltaScore'] = [maxscores[sn] - sc for sn,sc in zip(data['ScanNr'], data['Score'])]
charge_cols = [c for c in pred_cols if c.startswith('charge_')]
def countCharges(x): return np.sum(np.any(x.to_numpy(), axis = 0))
chargeCount = data.groupby('seq')[charge_cols].apply(countCharges).to_dict()
data['chargeCount'] = [chargeCount[s] for s in data['seq']]

data = data[['SpecId', 'Label', 'ScanNr']+pred_cols+['Peptide', 'Proteins']]

#run percolator on PSMs
data.to_csv('sipros.pin', sep = '\t', index = False)
# subprocess.run('percolator -U -m sipros.pout sipros.pin',
#                shell = True)
results = pd.read_csv('sipros.pout', sep = '\t')

#format results for further Sipros processing
results = results[results['q-value'] < 0.01]
results['scan'] = [i.split('@')[2] for i in results['PSMId']]
def best_per_scan(df):
    return df[df['posterior_error_prob'] == df['posterior_error_prob'].min()]
results = results.groupby('scan').apply(best_per_scan)
good_ids = set(results['PSMId'])
raw_data = raw_data[[i in good_ids for i in raw_data['SpecId']]]





