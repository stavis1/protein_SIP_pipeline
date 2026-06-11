#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:36:34 2026

@author: stavis
"""

import os
import re
import sys
import subprocess
import numpy as np
import pandas as pd
import configparser
from functools import cache
from collections import Counter

parser = configparser.ConfigParser()
parser.read(sys.argv[1])

#merge all individual label % searches into a single table
def read(f): return pd.read_csv(f, sep = '\t', comment = '#', dtype_backend="pyarrow")
raw_data = pd.concat([read(f) for f in os.listdir() if f.endswith('.sip')]).reset_index(drop = True)

#add columns needed for final output
raw_data['SpecId'] = ['@'.join(str(f) for f in fields) for fields in zip(raw_data['Filename'],
                                                                         raw_data['SearchName'],
                                                                         raw_data['ScanNumber'],
                                                                         raw_data['ParentCharge'],
                                                                         raw_data['Rank'])]
def deltaz(df):
    df = df.sort_values('Score')
    vals = list(df['Score'])
    dz = [vals[i]-vals[i-1] for i in range(1, len(vals))]
    return pd.Series([0] + dz, index = df.index)
raw_data['DeltaZ'] = raw_data.groupby('ScanNumber')[['Score']].apply(deltaz).droplevel('ScanNumber')
raw_data['DeltaP'] = [np.nan]*raw_data.shape[0]
raw_data['MassErrorDa'] = raw_data['CalculatedParentMass'] - raw_data['MeasuredParentMass']
raw_data['MassErrorPPM'] = raw_data['MassErrorDa'] * (raw_data['CalculatedParentMass']/1e6)
raw_data['ProteinCount'] = [len(p.split(',')) for p in raw_data['ProteinNames']]
raw_data['TargetMatch'] = ['F' if 'Rev_' in p else 'T' for p in raw_data['ProteinNames']]

#metadata columns needed for percolator
data = raw_data.copy()
data['Label'] = [1 if t == 'T' else -1 for t in data['TargetMatch']]
data['Peptide'] = data['IdentifiedPeptide']
data['Proteins'] = [('DECOY_' if 'Rev_' in p else '')+p for p in data['ProteinNames']]
data['ScanNr'] = data['ScanNumber']

#predictors for percolator
@cache
def getseq(p):
    return re.search(r'\[([^]]+)\]', p).group(1)
data['seq'] = [getseq(p) for p in data['IdentifiedPeptide']]
pred_cols = ['MeasuredParentMass',
             'Rank',
             'Score',
             'peplen',
             'N_missed',
             'MassErrorDa',
             'label_percent',
             'Nenz',
             'Cenz',
             'seqCount',
             'scanCount',
             'maxProtCount',
             'deltaScore',
             'DeltaZ',
             'chargeCount']
for charge in set(data['ParentCharge']):
    colname = f'charge_{charge}'
    pred_cols.append(colname)
    data[colname] = np.array(data['ParentCharge'] == charge, dtype = int)
data['peplen'] = [len(s) for s in data['seq']]
data['N_missed'] = [len(re.findall('[KR](?=.)', s)) for s in data['seq']]
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
init_tr = parser['Percolator']['train-fdr-initial']
tr = parser['Percolator']['trainFDR']
te = parser['Percolator']['testFDR']
subprocess.run(f'percolator -U --train-fdr-initial {init_tr} --trainFDR {tr} --testFDR {te} -m sipros.pout sipros.pin',
               shell = True)
results = pd.read_csv('sipros.pout', sep = '\t')

#filter raw results based on percolator output
results = results[results['q-value'] < 0.01]
results['scan'] = [i.split('@')[2] for i in results['PSMId']]
def best_per_scan(df):
    return df[df['posterior_error_prob'] == df['posterior_error_prob'].min()]
results = results.groupby('scan').apply(best_per_scan)
good_ids = set(results['PSMId'])
data = raw_data[[i in good_ids for i in raw_data['SpecId']]]

#format output tables
psm_cols = ['Filename',
            'ScanNumber',
            'ParentCharge',
            'MeasuredParentMass',
            'CalculatedParentMass',
            'MassErrorDa',
            'MassErrorPPM',
            'ScanType',
            'SearchName',
            'ScoringFunction',
            'Score',
            'DeltaZ',
            'DeltaP',
            'IdentifiedPeptide',
            'OriginalPeptide',
            'ProteinNames',
            'ProteinCount',
            'TargetMatch']
data = data[psm_cols]
filename = next(f for f in data['Filename'])
data.to_csv(f'{filename}.psm.txt', sep = '\t', index = False)

pep_cols = ['IdentifiedPeptide',
            'ParentCharge',
            'OriginalPeptide',
            'ProteinNames',
            'ProteinCount',
            'TargetMatch',
            'SpectralCount',
            'BestScore',
            'PSMs',
            'ScanType',
            'SearchName']
metadata = sorted(set(pep_cols).intersection(psm_cols), key = lambda c: pep_cols.index(c))
def make_peptide(df):
    first = next(i for i in df.index)
    row = df.loc[first, metadata]
    row['BestScore'] = max(df['Score'])
    row['PSMs'] = ','.join(f'{filename}[{s}]' for s in df['ScanNumber'])
    row['SpectralCount'] = df.shape[0]
    return row
peptides = data.groupby(['IdentifiedPeptide', 'ParentCharge'])[data.columns].apply(make_peptide)
peptides[pep_cols].to_csv(f'{filename}.pep.txt', sep = '\t', index = False)

