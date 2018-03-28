This supplement contains all of the scripts needed to reproduce the main results
of the manuscript.

To create the figures presented, the following simulations were ran using:

 - Intel Xeon 3.5 GHz processor with 8 cores and 32 gbs of RAM
 - The quadMINOS solver using the solvemepy (~4gb of RAM required)
 - Python 3.6

 Dependency versions
 - COBRApy v5.11.0
 - matplotlib v2.0.0
 - sympy v1.0.0
 - pandas v0.21.0
 - scipy v0.19.0
 - ecolime v0.0.7
 - cobrame v0.0.7
 - solvemepy v1.0.1


Note: Figure 4 and Table 3 are created using iJL1678b. If ecolime is installed, 
run build_me_model.py in the ecolime repository. This will used the information 
contained in this directory to reconstruct the iJL1678b model and output it as a
JSON. If ecolime is not installed, the model is already deposited in this 
supplementary directory. The remaining scripts in this supplemental directory 
will first search for the model built in ecolime before using the one already 
deposited.

Figure 5 is created using a version of iOL1650 created using COBRAme. It is
contained in this supplementary directory as COBRAme_iOL1650.json.


*************************** To create Figure 4 *******************************
1) If the quadMINOS solver and solvemepy is installed, run the run_qminos_simulation.py
file. This will solve the iJL1678b ME-model for each decimal point precision from
.1 to 1e-15 and track the time taken to solve. It will then use these growth rates
obtained at each precision and run flux variability analysis for the reactions
listed in the script as well as the transcription and translation reactions
required to produce the catalyzing enzyme. If the solver is not installed,
these files are deposited in already in the simulation_output folder. 

2) Run make_figure_4.py to plot the results from 1.

3) To create Figure S2, run make_figure_S2.py to plot the results
from 1 for two reactions in addition to PGI, which was used in figure 4.

*************************** To create Figure 5 ******************************
1) If the quadMINOS solver and solvemepy is installed, run the run_simulation.py
file. This will output the metabolic, transcription and translation fluxes
from iOL1650b into the COBRAme_simulations folder. If the solver is not
installed, the outputs are already deposited in the folder.

2) Run make_figure_5.py. This will run a pairwise comparison between the
deposited metabolic, transcription and translation fluxes in iOL1650_simulations
(these simulations are already deposited and were ran using the previous
ME-model) and those output from step 1. It will then output Figure 5.

*************************** To create Table 3 *******************************
1) If the quadMINOs solver and solvemepy is installed, run essentiality.py.
This will block the reactions synthesizing all 1678 genes in iJL1678b one-by-one.
It will then plug in .1 for mu and optimize for growth in glucose aerobic in
silico media. If the solution is infeasible then the gene is considered
essential and vice versa. This will output a file into the simulation_output
directory. If the solver is not installed the output is already deposited in
the output directory.

2) Run make_table_3.py. This will use a single gene knockout essentiality
data set (Monk_essentiality.csv from PMID: 29020004) to determine whether
the genes contained in both iOL1650 (based iOL1650_essentiality.xlsx, the
essentiality predictions from PMID: 24084808) and iJL1678b were correctly
predicted in each of the models. It will summarize these results in csv
spreadsheet.

************************** To create Figure S1 *****************************
1) If the quadMINOs solver and solvemepy is installed, run
run_qminos_simulations.py. This will sweep through maximum glucose uptake
rate values from -1 to -11 mmol/gDW/hr and optimize for growth rate with
ATPM as the objective. The simulations are output in the glucose_sweep directory.

2) Run make_figure_S1.py to plot the optimal growth rate for each
maximum uptake rate value in the top axis. The dual/shadow price values for 5
key metabolites/constraints are plotted for the same glucose uptake rates in
the bottom axis.
