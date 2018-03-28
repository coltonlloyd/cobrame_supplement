# coding: utf-8

from __future__ import print_function, division, absolute_import

import json
from os.path import abspath, dirname

# from solvemepy repository
from qminospy.me1 import ME_NLP1
from cobrame.io.json import load_json_me_model

here = dirname(abspath(__file__))


def solve_me_model(me, max_mu, precision=1e-14):

    # The object containing solveME methods--composite that uses a ME-model
    # object
    me_nlp = ME_NLP1(me, growth_key='mu')
    # Run bisection to find solution
    muopt, hs, xopt, cache = me_nlp.bisectmu(precision=precision, mumax=max_mu)
    me.solution.f = me.solution.x_dict['biomass_dilution']


def save_solution(model, solution=None):
    if not solution:
        solution = model.solution

    with open('%s/COBRAme_simulations/COBRAme_metabolic.json' % here,
              'w') as f:
        json.dump(model.get_metabolic_flux(solution=solution), f)

    with open('%s/COBRAme_simulations/COBRAme_transcription.json' % here,
              'w') as f:
        json.dump(model.get_transcription_flux(solution=solution), f)

    with open('%s/COBRAme_simulations/COBRAme_translation.json' % here,
              'w') as f:
        json.dump(model.get_translation_flux(solution=solution), f)


if __name__ == '__main__':
    # Load iOL1650 model constructed by COBRAme
    model = load_json_me_model('%s/COBRAme_iOL1650.json' % here)
    solve_me_model(model, 1., precision=1e-14)
    save_solution(model)
