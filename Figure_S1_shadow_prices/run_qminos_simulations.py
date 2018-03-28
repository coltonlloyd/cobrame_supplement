# coding: utf-8

from __future__ import print_function, division, absolute_import

from os.path import abspath, dirname
import json

import ecolime
from qminospy.me1 import ME_NLP1
import numpy as np
from cobrame.io.json import load_json_me_model

model_loc = '%s/me_models' % dirname(abspath(ecolime.__file__))

here = dirname(abspath(__file__))
saved_model_loc = dirname(here)
out_loc = '%s/glucose_sweep' % here

# Check for model in ecolime output directory first. If not present, load
# deposited model in supplement
try:
    me = load_json_me_model('%s/iJL1678b.json' % model_loc)
except:
    me = load_json_me_model('%s/iJL1678b.json' % saved_model_loc)


def run_glucose_sweep():
    objective = 'ATPM'
    uptake_values = np.linspace(1, 11, 20)
    me.objective = objective

    # Compile all mu expressions to speed up simulations
    me_nlp = ME_NLP1(me, growth_key='mu')
    me_nlp.compiled_expressions = me_nlp.compile_expressions()

    hs = None
    for uptake in uptake_values:
        print(uptake)
        me.reactions.EX_glc__D_e.lower_bound = -uptake

        # solve and reuse previous basis
        muopt, hs, xopt, cache = me_nlp.bisectmu(precision=1e-15,
                                                 mumax=1.5, basis=hs)
        me.solution.f = me.solution.x_dict['biomass_dilution']

        # Save with leading 0 if necessary
        out_file = '%06.3f_glucose_flux.json' % uptake
        with open('%s/%s' % (out_loc, out_file), 'w') as f:
            json.dump(me.solution.x_dict, f)

        out_dual_file = '%06.3f_glucose_dual.json' % uptake
        with open('%s/%s' % (out_loc, out_dual_file), 'w') as f:
            json.dump(me.solution.y_dict, f)


if __name__ == '__main__':
    run_glucose_sweep()
