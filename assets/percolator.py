#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 13:36:34 2026

@author: stavis
"""

import os
import re
# import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter#, defaultdict
# from statsmodels.nonparametric.smoothers_lowess import lowess

# predrt = pd.read_csv('predicted_rt.tsv')
# predscan = defaultdict(lambda:np.nan, {s:p for s,p in zip(predrt['seq'], predrt['predicted_rt'])})
# psms = pd.read_csv('SBF_boncat-facs-proteomics_13C_T1_R1_F1_consensus_fdr_filter.mzTab.psms',
#                    sep = '\t')
# psms['ScanNr'] = [int(re.search(r'scan=(\d+)', s).group(1)) for s in psms['spectra_ref']]
# psms['predrt'] = [predrt[s] for s in psms['sequence']]
# predscan = lowess(psms['ScanNr'], 
#                   psms['predrt'], 
#                   frac = 0.5,
#                   xvals = list(predrt.values()))
# predscan = {k:p for k,p in zip(predrt.keys(), predscan)}

def make_perc_table(file):
    data = pd.read_csv(file,
                       sep = '\t',
                       comment = '#')
    data['seq'] = [re.search(r'\[([^]]+)\]', p).group(1) for p in data['IdentifiedPeptide']]
    
    #metadata columns needed for percolator
    data['Label'] = [-1 if 'Rev_' in p else 1 for p in data['ProteinNames']]
    data['Peptide'] = data['IdentifiedPeptide']
    data['Proteins'] = [('DECOY_' if 'Rev_' in p else '')+p for p in data['ProteinNames']]
    data['ScanNr'] = data['ScanNumber']
    data['SpecId'] = ['@'.join(str(f) for f in fields) for fields in zip(data['Filename'],
                                                                         data['SearchName'],
                                                                         data['ScanNumber'],
                                                                         data['ParentCharge'],
                                                                         data['Rank'])]
    #predictors for percolator
    pred_cols = ['MeasuredParentMass',
                 'Rank',
                 'Score',
                 'peplen',
                 'N_missed',
                 'mass_err',
                 'label_percent',
                 # 'scan_err',
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
    data['label_percent'] = [int(re.search(r'_(\d+)Pct', file).group(1))/100000]*data.shape[0]
    # data['predscan'] = [predscan[s] for s in data['seq']]
    # data['scan_err'] = data['ScanNr'] - data['predscan']
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
    return data

#run percolator on PSMs
sips = [f for f in os.listdir() if f.endswith('.sip')]
data = pd.concat([make_perc_table(f) for f in sips])
data.to_csv('sipros.pin', sep = '\t', index = False)
# subprocess.run('percolator -U -m sipros.pout sipros.pin',
#                shell = True)

# =============================================================================
# testing
# =============================================================================
if False:    
    targets = data[data['Label'] == 1]
    decoys = data[data['Label'] == -1]
    for col in data.columns[3:-2]:
        fig, ax = plt.subplots(dpi = 900, figsize = (5,5))
        
        if len(set(data[col])) < 5:
            t_counts = sorted(Counter(targets[col]).items())
            d_counts = sorted(Counter(decoys[col]).items())
            ax.bar([c[0] - 0.15 for c in t_counts],
                   [c[1] for c in t_counts],
                   color = 'g',
                   width = 0.3,
                   label = 'Target',
                   alpha = 0.5)
            ax.bar([c[0] + 0.15 for c in d_counts],
                   [c[1] for c in d_counts],
                   color = 'k',
                   width = 0.3,
                   label = 'Decoy',
                   alpha = 0.5)
        else:
            bins = np.linspace(min(data[col]),
                               max(data[col]),
                               200)
            ax.hist(targets[col], 
                    bins = bins,
                    color = 'g',
                    label = 'Target',
                    alpha = 0.5)
            ax.hist(decoys[col], 
                    bins = bins,
                    color = 'k',
                    label = 'Decoy',
                    alpha = 0.5)
        ax.legend()
        ax.set_title(col)

    # def mods(text):
    #     locs = [i+1 for i,c in enumerate(text) if c == 'C']
    #     return '|'.join(str(e) for m  in zip(locs, ['Carbamidomethyl']*len(locs)) for e in m)
    
    seqs = list({re.search(r'\[([^]]+)\]', p).group(1) for p in data['Peptide']})
    deeplc = pd.DataFrame({'seq':seqs,
                           'modifications':['']*len(seqs),
                           'tr':['']*len(seqs)})
    deeplc.to_csv('peptides.csv', index = False)

    # fig, ax = plt.subplots(dpi = 900, figsize = (5,5))
    # ax.scatter(psms['ScanNr'], [predscan[s] for s in psms['sequence']], s = 1, c = 'k', marker = '.')
    
    # cal_data = psms.groupby('sequence')['ScanNr'].apply(np.mean)
    # deeplc_calibrate = pd.DataFrame({'seq':cal_data.index,
    #                                  'modifications':['']*len(cal_data),
    #                                  'tr':list(cal_data)})
    # deeplc_calibrate.to_csv('deeplc_calibrate.csv', index = False)
    
    