# Code extracted from Jupyter Notebook cell

# import libraries 
import mineapy
from cobra.io import load_matlab_model,load_json_model
import pandas as pd
from mineapy.core.taskEnrich import TaskEnrichment
from mineapy.core.thermo_model import ThermoModel_WithoutInfo
from mineapy.core.rxnExp import ReactionExp
import sys
import os
import pickle
# Add the path to the tutorials directory to sys.path
sys.path.append(os.path.abspath(os.path.join('..')))

#print(os.path.abspath(os.path.join('..')))

from data_utility import getLiverStage



# Code extracted from Jupyter Notebook cell

# This function load the gene expression data
df=getLiverStage()
print(df.head())



# Code extracted from Jupyter Notebook cell

# load the metabolic model and Python version is saved in the models folder. 
ipbe_blood=pickle.load(open('../models/ipbeblood_py.pickle','rb'))
setattr(ipbe_blood, 'annotation',ipbe_blood._annotation)
for item in ipbe_blood.metabolites:
    setattr(item, 'annotation',item._annotation)
for item in ipbe_blood.reactions:
    setattr(item, 'annotation',item._annotation)
for item in ipbe_blood.genes:
    setattr(item, 'annotation',item._annotation)


## get biomass reaction 
biomass = ipbe_blood.reactions.get_by_id('biomass')
biomass.lower_bound=0
tfa_solution = ipbe_blood.optimize()
tfa_value = tfa_solution.objective_value
print('Solution found : {0:.5g}'.format(tfa_value))


# print genes
gene_ids=[g.id for g in ipbe_blood.genes]
print(gene_ids)



# Code extracted from Jupyter Notebook cell

# A value of high_cutoff=0.15 indicates that the top 15% of highly expressed genes.
# and the bottom 15% of lowly expressed genes will be selected when low_cutoff=0.15.
gene_exp={'gene_id':df['Gene ID'].to_list(),'exp_val':df['Mean_RPKM'].to_list(),'high_cutoff':0.15,'low_cutoff':0.15}

exp_analysis=ReactionExp(ipbe_blood,gene_exp=gene_exp)

#params_rxns={'up_rxns':reg_analysis.up_rxns,'down_rxns':reg_analysis.down_rxns}
params_rxns={'high_rxns':exp_analysis.high_rxns,'low_rxns':exp_analysis.low_rxns}




# Code extracted from Jupyter Notebook cell

sol=ipbe_blood.optimize()
print(sol,type(ipbe_blood.solver))



# Code extracted from Jupyter Notebook cell

path_to_params = path_to_params = '../input/Minea_parameter_ipbe.yaml'
task_enrich = TaskEnrichment(ipbe_blood,path_to_params,params_rxns)
task_enrich.run()



# Code extracted from Jupyter Notebook cell



