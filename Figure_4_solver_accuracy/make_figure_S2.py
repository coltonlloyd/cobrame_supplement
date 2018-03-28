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


def make_figure(fva, ax):

    # create array with precision of each FVA simulation
    prec = []
    for info in fva:
        prec.append(-log10(info[0]))

    # create arrays of max in min values of each reaction from FVA output
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
            label = rxn.split('_')[0].replace(plot_rxn, 'metabolic')
            if 'REV' in rxn:
                label += ' (Reverse) %i' % metabolic_rev_num
                metabolic_rev_num += 1
                marker = '^'
            else:
                label += ' %i' % metabolic_num
                metabolic_num += 1
                marker = None

        label = label.title()  # capitalize first letter

        # plot both min and max values for each reaction
        line = ax.plot(prec, values, label=label, marker=marker)
        ax.plot(prec, out_dict[rxn + '_min'], '--',
                color=line[-1].get_color(), marker=marker)

    # Format table and add labels
    ax.locator_params('x', nbins=10)
    ax.set_ylim([1e-16, 100])
    ax.set_xlim([0, 16])
    ax.set_ylabel(r'Max and Min Flux ($\frac{mmol}{gDW \cdot hr}$)',
                  fontsize=20)
    fig.text(.5, .025, 'Solve Precision (Number of Decimal Points)',
             fontsize=20, ha='center')
    ax.set_title('Flux Variability \n of %s Reactions' % plot_rxn, fontsize=20)


if __name__ == '__main__':
    here = dirname(abspath(__file__))
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharex=True, sharey=True)
    fig.subplots_adjust(bottom=0.15, wspace=.3)

    # Make two plots similar to the main figure with MGSA and E4PD reactions
    for i, plot_rxn in enumerate(['MGSA', 'E4PD']):

        out_loc = '%s/simulation_output' % here

        with open('%s/%s_fva.json' % (out_loc, plot_rxn), 'r') as f:
            fva = json.load(f)

        make_figure(fva, axes[i])

        # Add and format legend for each subplot
        lgd = axes[i].legend(fontsize=14, ncol=2, bbox_to_anchor=(.5, -.2),
                             loc='upper center',
                             title='%s Reactions' % plot_rxn)
        lgd.get_title().set_fontsize('17')

    fig.savefig('%s/Figure_S2.png' % here,
                bbox_extra_artists=(lgd,), bbox_inches='tight')
