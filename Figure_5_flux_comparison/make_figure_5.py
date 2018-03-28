# coding: utf-8

from __future__ import print_function, division, absolute_import

import json
from os.path import abspath, dirname
from cycler import cycler

import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 20
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
plt.rcParams['legend.fontsize'] = 'xx-large'
plt.rcParams['axes.prop_cycle'] = cycler(color=['#0099E6', '#F23814',
                                                '#DFB78B', '#B6C3C5',
                                                '#12255A'])
plt.rcParams.update({'figure.autolayout': True})

textsize = 25
markersize = 80

here = dirname(abspath(__file__))

iol_sims_loc = '%s/iOL1650_simulations' % here
cobrame_sims_loc = '%s/COBRAme_simulations' % here


if __name__ == '__main__':

    # ****************** Get translation fluxes *******************************
    # Get translation fluxes and create dataframe for comparison
    with open('%s/COBRAme_translation.json' % cobrame_sims_loc, 'r') as f:
        cobrame_translation_flux = json.load(f)
    with open('%s/iOL1650_translation.json' % iol_sims_loc, 'r') as f:
        iol1650_translation_flux = json.load(f)

    translation_compare_dict = {'COBRAme': cobrame_translation_flux,
                                'iOL1650': iol1650_translation_flux}
    translation_compare_df =\
        pd.DataFrame.from_dict(translation_compare_dict).dropna(how='any')

    # ****************** Get transcription fluxes *****************************
    with open('%s/COBRAme_transcription.json' % cobrame_sims_loc, 'r') as f:
        cobrame_transcription_flux = json.load(f)
    with open('%s/iOL1650_transcription.json' % iol_sims_loc, 'r') as f:
        iol1650_transcription_flux = json.load(f)

    transcription_compare_dict = {'COBRAme': cobrame_transcription_flux,
                                  'iOL1650': iol1650_transcription_flux}
    transcription_compare_df =\
        pd.DataFrame.from_dict(transcription_compare_dict).dropna(how='any')

    # ****************** Get metabolic reaction fluxes ************************
    with open('%s/COBRAme_metabolic.json' % cobrame_sims_loc, 'r') as f:
        cobrame_metabolic_flux = json.load(f)
    with open('%s/iOL1650_metabolic.json' % iol_sims_loc, 'r') as f:
        iol1650_metabolic_flux = json.load(f)

    metabolic_compare_dict = {'COBRAme': cobrame_metabolic_flux,
                              'iOL1650': iol1650_metabolic_flux}
    metabolic_compare_df =\
        pd.DataFrame.from_dict(metabolic_compare_dict).dropna(how='any')

    # Copy dataframe with negative values to plot on linear axis
    met_df = metabolic_compare_df.copy()

    # Set values < 1e-14 to 1e-14 (This will only include ractions with 0 flux)
    transcription_compare_df = transcription_compare_df.clip_lower(1e-14)
    translation_compare_df = translation_compare_df.clip_lower(1e-14)
    metabolic_compare_df = abs(metabolic_compare_df).clip_lower(1e-14)

    # Get correlations of fluxes predicted by iOL1650 vs COBRAme model
    transc_corr = transcription_compare_df.corr().values[1][0] ** 2
    transl_corr = translation_compare_df.corr().values[1][0] ** 2
    met_corr = metabolic_compare_df.corr().values[1][0] ** 2

    # ******* Plot all fluxes on one plot using log log scale ***************
    fig = plt.figure(figsize=(15, 10))
    ax = plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=3)
    ax.figure.subplots_adjust(bottom=0.15, wspace=.3)
    ax.loglog()

    ax.scatter(transcription_compare_df['COBRAme'],
               transcription_compare_df['iOL1650'], s=markersize, c='#0099E6',
               label=r'Transcription Flux ($R^2$ = %.3f)' % transc_corr)

    ax.scatter(translation_compare_df['COBRAme'],
               translation_compare_df['iOL1650'], s=markersize, c='#DFB78B',
               label='Translation Flux ($R^2$ = %.3f)' % transl_corr)

    ax.scatter(metabolic_compare_df['COBRAme'],
               metabolic_compare_df['iOL1650'], s=markersize, c='#12255A',
               label='Metabolic Flux ($R^2$ = %.3f)' % met_corr)

    # ************** Format axes and add labels/legend ************************
    ax.set_ylabel(r'$i$OL1650 ($\frac{mmol}{gDW \cdot hr}$)',
                  fontsize=textsize, color='k')
    fig.text(.5, -.05,
             r'$i$OL1650, '
             r'constructed with COBRAme ($\frac{mmol}{gDW \cdot hr}$)',
             fontsize=textsize, ha='center')
    ax.set_xlim([1e-14, 1000])
    ax.set_ylim([1e-14, 1000])
    ax.plot([1e-14, 1000], [1e-14, 1000], c='#F23814')
    ax.legend(loc='lower right', fontsize=20, markerscale=1.5,
              handletextpad=-.1)
    ax.set_title('All ME-model Reaction Fluxes', fontsize=textsize)

    ax.axes.xaxis.set_tick_params(labelcolor='k')
    ax.axes.yaxis.set_tick_params(labelcolor='k')

    # ************ Plot predicted fluxes on separate linear axes *************
    # Plot transcription fluxes
    ax0 = plt.subplot2grid((3, 3), (0, 2))
    ax0 = transcription_compare_df.plot(kind='scatter', x='COBRAme',
                                        y='iOL1650',
                                        loglog=False, xlim=[0, 6e-6],
                                        ylim=[0, 6e-6], ax=ax0, s=markersize,
                                        c='#0099E6')
    ax0.plot([0, .1], [0, .1], c='#F23814')

    # Format axis
    ax0.set_xlabel('')
    ax0.set_ylabel('')
    ax0.set_title('Transcription', fontsize=textsize)
    ax0.axes.xaxis.set_tick_params(labelcolor='k')
    ax0.axes.yaxis.set_tick_params(labelcolor='k')
    ax0.locator_params(nbins=3)
    ax0.ticklabel_format(style='sci', scilimits=(0, 0))
    ax0.get_yaxis().get_offset_text().set_position((-.4, .1))

    # Plot transcription fluxes
    ax1 = plt.subplot2grid((3, 3), (1, 2))
    ax1 = translation_compare_df.plot(kind='scatter', x='COBRAme', y='iOL1650',
                                      loglog=False, xlim=[0, 6e-4],
                                      ylim=[0, 6e-4], ax=ax1, s=markersize,
                                      c='#DFB78B')
    ax1.plot([0, .1], [0, .1], c='#F23814')

    # Format axis
    ax1.minorticks_on()
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    ax1.set_title('Translation', fontsize=textsize)
    ax1.axes.xaxis.set_tick_params(labelcolor='k')
    ax1.axes.yaxis.set_tick_params(labelcolor='k')
    ax1.locator_params(nbins=3)
    ax1.ticklabel_format(style='sci', scilimits=(0, 0))
    ax1.get_yaxis().get_offset_text().set_position((-.4, .1))

    # Plot metabolic fluxes
    ax2 = plt.subplot2grid((3, 3), (2, 2))
    ax2 = met_df.plot(kind='scatter', x='COBRAme', y='iOL1650',
                      ylim=[-25, 25], xlim=[-25, 25], ax=ax2, s=markersize,
                      c='#12255A')
    ax2.plot([-60, 60], [-60, 60], c='#F23814')

    # Format axis
    ax2.set_xlabel('')
    ax2.set_ylabel('')
    ax2.set_title('Metabolic', fontsize=textsize)
    ax2.axes.xaxis.set_tick_params(labelcolor='k')
    ax2.axes.yaxis.set_tick_params(labelcolor='k')
    ax2.locator_params(nbins=2)

    # Save figure
    ax.figure.savefig('%s/Figure_5.png' % here, bbox_inches="tight")
