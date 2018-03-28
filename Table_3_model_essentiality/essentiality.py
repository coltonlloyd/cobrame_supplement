from __future__ import print_function, division, absolute_import

import re

import json
from os.path import abspath, dirname

import ecolime
from cobrame.core.reaction import *
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


def compute_gene_essentiality_at_gr(gr):
    me_nlp = ME_NLP1(me, growth_key='mu')
    expressions = me_nlp.compile_expressions()
    me_nlp.compiled_expressions = expressions

    hs = None

    all_genes = me.metabolites.query(re.compile("^RNA_b[0-9]"))

    results = {}
    for gene_RNA in list(all_genes):

        default_bounds = {}
        for r in gene_RNA.reactions:
            if not r.id.startswith("DM") and not \
                    isinstance(r, TranscriptionReaction):
                default_bounds[r] = (r.lower_bound, r.upper_bound)
                r.knock_out()
        x, status, hs = me_nlp.solvelp(gr, basis=hs)

        if status == 'optimal':
            results[gene_RNA.id] = 'NONESSENTIAL'
        else:
            results[gene_RNA.id] = 'ESSENTIAL'

        print("%s\t%s" % (gene_RNA.id.split("_")[1], str(status)))

        with open("%s/iJL1678b_essentiality_%.2f_gr.json" % (out_loc, gr),
                  "w") as outfile:
            json.dump(results, outfile, indent=True)

        # Reset bounds
        for r in default_bounds:
            r.lower_bound = default_bounds[r][0]
            r.upper_bound = default_bounds[r][1]

if __name__ == '__main__':
    compute_gene_essentiality_at_gr(.1)
