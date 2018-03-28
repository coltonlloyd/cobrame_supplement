from __future__ import print_function, division, absolute_import

from cycler import cycler
from glob import glob
import pickle
from os.path import dirname, abspath

from matplotlib import pyplot as plt
import pandas as pd

plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 25
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['grid.color'] = 'white'
plt.rcParams['figure.edgecolor'] = 'black'
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['axes.grid'] = False
plt.rcParams['legend.edgecolor'] = 'k'
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.markersize'] = 8

plt.rcParams['legend.fontsize'] = 'large'
plt.rcParams['axes.prop_cycle'] = cycler(color=['#0099E6', '#F23814',
                                                '#DFB78B', '#B6C3C5',
                                                '#12255A'])

here = dirname(abspath(__file__))

if __name__ == '__main__':
    y = []
    x = []

    # Create DataFrame containing shadow prices of all metabolites for each
    # glucose uptake rate in glucose sweep
    df = pd.DataFrame()
    for i in sorted(glob('%s/glucose_sweep/*' % here)):
        with open(i, 'rb') as f:
            sol = pickle.load(f)
        y.append(sol.f)
        x.append(float(i.split('/')[-1].split('_')[0]))
        df_new = pd.DataFrame.from_dict(sol.y_dict, orient='index')
        df_new.columns = [i.split('/')[-1].split('_')[0]]
        df = df.join(df_new, how='outer')

    fig, axes = plt.subplots(2, 1, figsize=(6., 7.5), sharex='col')

    axes[0].plot(x, y, linewidth=3)
    axes[0].set_ylabel(r'Growth rate ($\mu^{-1}$)')
    mets = ['protein_biomass', 'biomass', 'mRNA_biomass', 'ribosome',
            'protein_dummy']
    for met in mets:
        axes[1].plot([float(i) for i in df.columns], df.T[met].values,
                     label=met, linewidth=3)

    axes[1].set_yscale('symlog')
    axes[1].set_xlabel(r'Glucose Uptake Rate($\frac{mmol}{gDW \cdot hr}$)')
    axes[1].set_ylabel('Shadow Price ' +
                       r'($\mathrm{\frac{mol_{ATP}}{mol_{metabolite}}}$)')
    axes[1].set_yticks([-1e6, -1e3, 0.0, 1e3, 1e6])

    plt.legend(loc=(.02, .35))
    fig.tight_layout()
    fig.savefig('%s/Figure_S1.png' % here)
