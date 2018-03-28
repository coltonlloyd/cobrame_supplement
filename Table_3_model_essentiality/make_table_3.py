
# coding: utf-8
from __future__ import print_function, division, absolute_import

from os.path import abspath, dirname
import math

import pandas as pd
import json

here = dirname(abspath(__file__))

# Load computational essentiality datasets
iOL1650_essentiality_df = pd.read_excel('%s/iOL1650_essentiality.xlsx' % here,
                                        index_col=0)
with open('%s/simulation_output/iJL1678b_essentiality_0.10_gr.json' % here,
          'r') as f:
    iJL1678b_essential_json = json.load(f)

# Load experimental essentiality dataset
df_essentiality_monk = pd.read_csv('%s/Monk_essentiality.csv' % here,
                                   index_col=1)

out = dict(TP=0, FP=0, FN=0, TN=0)
out_iol = dict(TP=0, FP=0, FN=0, TN=0)

for gene in iOL1650_essentiality_df.index:

    # These are LB essential genes
    if gene not in df_essentiality_monk.index:
        essentiality = 'ESSENTIAL'
    else:
        essentiality = \
            'ESSENTIAL' if (float(df_essentiality_monk.loc[gene, 'growth'])
                            < .5) else 'NONESSENTIAL'
    if 'RNA_' + gene in iJL1678b_essential_json:
        iML1678b_essentiality = iJL1678b_essential_json['RNA_' + gene]
    else:
        continue

    # iJL1678b
    if essentiality == 'ESSENTIAL' and iML1678b_essentiality == 'ESSENTIAL':
        out['TN'] += 1
    elif essentiality == 'ESSENTIAL' and iML1678b_essentiality != 'ESSENTIAL':
        out['FP'] += 1
    elif essentiality == 'NONESSENTIAL' and \
            iML1678b_essentiality == 'NONESSENTIAL':
        out['TP'] += 1
    elif essentiality == 'NONESSENTIAL' and \
            iML1678b_essentiality != 'NONESSENTIAL':
        out['FN'] += 1

    # iOL1650
    iOL1650_essentiality = iOL1650_essentiality_df.loc[gene, 'ME essentiality']
    if essentiality == 'ESSENTIAL' and iOL1650_essentiality == 'ESSENTIAL':
        out_iol['TN'] += 1
    elif essentiality == 'ESSENTIAL' and iOL1650_essentiality != 'ESSENTIAL':
        out_iol['FP'] += 1
    elif essentiality == 'NONESSENTIAL' and \
            iOL1650_essentiality == 'NONESSENTIAL':
        out_iol['TP'] += 1
    elif essentiality == 'NONESSENTIAL' and \
            iOL1650_essentiality != 'NONESSENTIAL':
        out_iol['FN'] += 1

out_df = pd.DataFrame([out, out_iol],
                      index=['iJL1678b_essentiality', 'iOL1650_essentiality'])

for i in out_df.index:
    TP = out_df.loc[i, 'TP']
    TN = out_df.loc[i, 'TN']
    FP = out_df.loc[i, 'FP']
    FN = out_df.loc[i, 'FN']
    mcc = ((TP*TN)-(FP*FN))/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    out_df.loc[i, 'Matthews Correlation Coefficient'] = mcc

out_df.to_csv('%s/Table_3.csv' % here)
