# coding: utf-8

from __future__ import print_function, division, absolute_import

import time
import json
import os.path
from os.path import abspath, dirname

import ecolime
from qminospy.me1 import ME_NLP1
from cobrame.io.json import load_json_me_model

model_loc = '%s/me_models' % dirname(abspath(ecolime.__file__))

here = dirname(abspath(__file__))
saved_model_loc = dirname(here)
out_loc = '%s/simulation_output' % here

# Check for model in ecolime output directory first. If not present, load
# deposited model in supplement
try:
    me = load_json_me_model('%s/iJL1678b.json' % model_loc)
except:
    me = load_json_me_model('%s/iJL1678b.json' % saved_model_loc)


def add_transcription_rxns(bnum, sim_reaction_set):
    for r in me.metabolites.get_by_id('RNA_'+bnum).reactions:
        if me.metabolites.get_by_id('RNA_'+bnum) in r.products:
            sim_reaction_set.add(r.id)


def add_all_rxns(rxn, sim_reaction_set):
    d = me.process_data.get_by_id(rxn)
    for r in d.parent_reactions:
        sim_reaction_set.add(r.id)
        if 'SPONT' in r.id:
            cp = 'CPLX_dummy'
        else:
            cp = r.id.split('_FWD_')[-1].split('_REV_')[-1]
        for tr in me.process_data.get_by_id(cp).stoichiometry:
            protein = tr.replace('protein_', '')
            bnum = protein.split('_')[0]
            sim_reaction_set.add('translation_' + bnum)
            add_transcription_rxns(bnum, sim_reaction_set)


def run_precision_solves(me, rxn):
    times_list = []
    for precision in [.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9,
                      1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15]:

        sim_reaction_set = set()
        add_all_rxns(rxn, sim_reaction_set)

        start = time.time()
        me_nlp = ME_NLP1(me, growth_key='mu')
        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision,
                                                 mumax=1.5)
        me.solution.f = me.solution.x_dict['biomass_dilution']

        end = time.time()
        total = end-start
        gr = me.solution.f
        a = (precision, total, gr)
        times_list.append(a)

        with open('%s/%s_time.json' % (out_loc, rxn), 'w') as f:
            json.dump(times_list, f)

    return times_list


def run_fva(me, rxn, times_list):
    fva_list = []
    hs = None
    for entry in times_list:

        precision = entry[0]
        gr = entry[2]

        me_nlp = ME_NLP1(me, growth_key='mu')

        sim_reaction_set = set()

        add_all_rxns(rxn, sim_reaction_set)
        print(sim_reaction_set)
        out, _ = me_nlp.varyme(gr, list(sim_reaction_set), basis=hs)
        hs = me_nlp.hs

        fva_list.append((precision, out, gr))
        with open('%s/%s_fva.json' % (out_loc, rxn), 'w') as f:
            json.dump(fva_list, f)


if __name__ == '__main__':
    # Run FVA for the following reactions
    rxns = ['dummy_reaction', 'AKGt2rpp', 'E4PD', 'PGI', 'TYRL', 'PSERT',
            'DHAD1', 'CYSTL', 'ENO', 'MGSA']
    for rxn in rxns:
        print(rxn)

        times_output = '%s/%s_time.json' % (out_loc, rxn)

        if os.path.isfile(times_output):
            print('using existing simulations')
            with open(times_output, 'r') as f:
                times = json.load(f)
        else:
            times = run_precision_solves(me, rxn)

        run_fva(me, rxn, times)
