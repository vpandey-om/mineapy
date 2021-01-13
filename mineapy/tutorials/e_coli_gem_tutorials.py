import mineapy
from cobra.io import load_matlab_model,load_json_model
import pandas as pd
from mineapy.core.taskEnrich import TaskEnrichment
from mineapy.core.thermo_model import ThermoModel_WithoutInfo
from mineapy.core.rxnExp import ReactionExp
import pytfa
from pytfa.io import  read_compartment_data, apply_compartment_data, \
    read_lexicon, annotate_from_lexicon
from pytfa.io.base import load_thermoDB




from pytfa.thermo.tmodel import ThermoModel


from os.path import join


cobra_model = load_json_model('./models/iJO1366.json')
genes=[g.id for g in cobra_model.genes]



path_to_params = './input/task_enrichment_gem_params.yaml'

context_df=pd.read_csv('./input/context.txt',sep='\t')
condition_df=pd.read_csv('./input/condition.txt',sep='\t')
### working with condition comparisons
gene_reg={'gene_id':condition_df['geneid'].to_list(),'fold_change':condition_df['fold change'].to_list(),'up_cutoff':1.35,'down_cutoff':float(1/2.5)}

reg_analysis=ReactionExp(cobra_model,gene_reg=gene_reg)

gene_exp={'gene_id':context_df['geneid'].to_list(),'exp_val':context_df['exp_val'].to_list(),'high_cutoff':0.15,'low_cutoff':0.15}

exp_analysis=ReactionExp(cobra_model,gene_exp=gene_exp)

#params_rxns={'up_rxns':reg_analysis.up_rxns,'down_rxns':reg_analysis.down_rxns}
params_rxns={'high_rxns':exp_analysis.high_rxns,'low_rxns':exp_analysis.low_rxns}

# 1)  enrichment based on cobra model
#import pdb; pdb.set_trace()
#cobra_model.solver='optlang-gurobi'
#cobra_model.solver= 'optlang-cplex'

# task_enrich = TaskEnrichment(cobra_model,path_to_params,params_rxns)
#
# task_enrich.run()


# 2) enrichment based on cobra model


# Paths



# Paths
path_to_model = join('.','models','iJO1366.json')
thermoDB = join('.','input','thermo_data.thermodb')
path_to_lexicon = join('.','models','iJO1366','lexicon.csv')
path_to_compartment_data = join('.','models','iJO1366','compartment_data.json')

# FBA
model = load_json_model(path_to_model)

fba_solution = model.optimize()
fba_value = fba_solution.objective_value

# Thermo prep

thermo_data = load_thermoDB(thermoDB)
lexicon = read_lexicon(path_to_lexicon)
compartment_data = read_compartment_data(path_to_compartment_data)

# Initialize the cobra_model
tfa_model = ThermoModel(thermo_data, model)
tfa_model.name = 'Tutorial'

# Annotate the cobra_model
annotate_from_lexicon(tfa_model, lexicon)
apply_compartment_data(tfa_model, compartment_data)

tfa_model.prepare()
tfa_model.convert()

# tfa_model.solver.configuration.verbosity = True
tfa_model.logger.setLevel = 30

tfa_solution = tfa_model.optimize()
tfa_value = tfa_solution.objective_value

# It might happen that the model is infeasible. In this case, we can relax
# thermodynamics constraints:

if tfa_value < 0.1:
    from pytfa.optim.relaxation import relax_dgo

    biomass_rxn = 'Ec_biomass_iJO1366_WT_53p95M'
    tfa_model.reactions.get_by_id(biomass_rxn).lower_bound = 0.9 * fba_value
    relaxed_model, slack_model, relax_table = relax_dgo(tfa_model, in_place=True)

    original_model, tfa_model = tfa_model, relaxed_model

    print('Relaxation: ')
    print(relax_table)

    tfa_solution = tfa_model.optimize()
    tfa_value = tfa_solution.objective_value


# Thermo prep


# biomass_rxn = 'Ec_biomass_iJO1366_WT_53p95M'
tfa_model.reactions.get_by_id(biomass_rxn).lower_bound = 0
tfa_model.solver= 'optlang-cplex'

task_enrich = TaskEnrichment(tfa_model,path_to_params,params_rxns)

task_enrich.run()
