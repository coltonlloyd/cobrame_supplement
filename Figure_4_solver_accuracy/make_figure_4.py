from __future__ import print_function, division, absolute_import

import glob
import json
from os.path import abspath, dirname
from math import log10
from collections import defaultdict

from matplotlib import pyplot as plt
from cycler import cycler
import pandas as pd

plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
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


def make_figure(fva):
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True)
    fig.subplots_adjust(bottom=0.15, wspace=.3)

    # get an array of the precisions computed
    # and create a dataframe w/ the time for each precision for each reaction
    prec = []
    comp_time_df = pd.DataFrame()
    for i, sim in enumerate(glob.glob('%s/*_time.json' % (out_loc))):
        with open(sim, 'r') as f:
            times = json.load(f)
        for info in reversed(times):
            if i == 0:
                prec.append(-log10(info[0]))
            comp_time_df.loc[i, -log10(info[0])] = info[1]

    # plot the times w/ standard deviation errorbars for each precision
    # solvetime
    ax = axes[1]
    ax.errorbar(comp_time_df.columns, comp_time_df.mean(), fmt='o',
                yerr=comp_time_df.std(), ecolor='k', capsize=10, elinewidth=3)
    ax.set_ylabel('Solve Time (Seconds)', fontsize=20)
    ax.set_title('Solve Time vs. Solver Precision', fontsize=20)

    # create arrays of max in min values of each reaction from FVA output
    rev_prec = list(reversed(list(prec)))
    out_dict = defaultdict(list)
    for i, entry in enumerate(fva):
        for key, value in entry[1].items():

            # 1e-15 is close to the lowest value floating double precision
            # can handle
            if value['minimum'] < 1e-15:
                min_value = 1e-15
            else:
                min_value = value['minimum']

            if value['maximum'] < 1e-15:
                max_value = 1e-15
            else:
                max_value = value['maximum']
            out_dict[key + '_min'].append(min_value)
            out_dict[key].append(max_value)

    # Loop through FVA reactions and plot
    ax = axes[0]
    transcription_num = translation_num = metabolic_num = metabolic_rev_num = 1
    already_ploted = []
    for rxn in sorted(out_dict):
        values = out_dict[rxn]
        if rxn in already_ploted:
            continue
        if '_min' in rxn:
            continue
        ax.semilogy()

        # Get label and marker type for each reaction category
        if 'transcription_' in rxn:
            label = 'Transcription %i' % transcription_num
            transcription_num += 1
            marker = 'x'
        elif 'translation_' in rxn:
            label = 'Translation %i' % translation_num
            translation_num += 1
            marker = 'o'
        else:
            label = rxn.split('_')[0].replace(plot_rxn, 'Metabolic')
            if 'REV' in rxn:
                label += ' (Reverse) %i' % metabolic_rev_num
                metabolic_rev_num += 1
                marker = '^'
            else:
                label += ' %i' % metabolic_num
                metabolic_num += 1
                marker = None

        # plot both min and max values for each reaction
        line = ax.plot(rev_prec, values, label=label, marker=marker)
        ax.plot(rev_prec, out_dict[rxn + '_min'], '--',
                color=line[-1].get_color(), marker=marker)

    # Format table and add labels / legend
    ax.locator_params('x', nbins=10)
    ax.set_ylim([1e-16, 100])
    ax.set_xlim([0, 16])
    ax.set_ylabel(r'Max and Min Flux ($\frac{mmol}{gDW \cdot hr}$)',
                  fontsize=20)
    fig.text(.5, .025, 'Solve Precision (Number of Decimal Points)',
             fontsize=20, ha='center')
    ax.set_title('Flux Variability \n of %s Reactions' % plot_rxn, fontsize=20)
    lgd = ax.legend(fontsize=14, ncol=2, bbox_to_anchor=(1, -.2),
                    loc='upper center', title='%s Reactions' % plot_rxn)
    lgd.get_title().set_fontsize('17')

    fig.savefig('%s/Figure_4_%s.png' % (here, plot_rxn),
                bbox_extra_artists=(lgd,), bbox_inches='tight')


if __name__ == '__main__':
    here = dirname(abspath(__file__))
    out_loc = '%s/simulation_output' % here

    for plot_rxn in ['PGI']:

        with open('%s/%s_time.json' % (out_loc, plot_rxn), 'r') as f:
            times = json.load(f)

        with open('%s/%s_fva.json' % (out_loc, plot_rxn), 'r') as f:
            fva = json.load(f)

        make_figure(fva)
