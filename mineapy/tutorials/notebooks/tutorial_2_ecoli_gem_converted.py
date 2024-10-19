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
from data_utility import getErythromycin_data



# Code extracted from Jupyter Notebook cell

## load the gene expression data for gene expression   Erythromycin 
ery_df=getErythromycin_data()

ery_df['significant(ERYvsETH)'].unique()
# Convert the 'log2FoldChange(ERYvsETH)' column to numeric, coercing errors to NaN
ery_df['log2FoldChange(ERYvsETH)'] = pd.to_numeric(ery_df['log2FoldChange(ERYvsETH)'], errors='coerce')

## find up and down regulated genes 
up_down_df=ery_df[(ery_df['significant(ERYvsETH)']=='UP')|(ery_df['significant(ERYvsETH)']=='DOWN')]
up_down_df['fold_change']=2 ** up_down_df['log2FoldChange(ERYvsETH)']

print(up_down_df)



# Code extracted from Jupyter Notebook cell

# load core metabolic genome_scale_model
cobra_model = load_json_model('../models/iJO1366.json')
genes=[g.id for g in cobra_model.genes]
sol=cobra_model.optimize()
print(sol,type(cobra_model.solver))



# Code extracted from Jupyter Notebook cell

### working with condition comparisons
gene_reg={'gene_id':up_down_df['Gene_id'].to_list(),'fold_change':up_down_df['fold_change'].to_list(),'up_cutoff':1.1,'down_cutoff':0.9}

reg_analysis=ReactionExp(cobra_model,gene_reg=gene_reg)
reg_analysis
params_rxns={'up_rxns':reg_analysis.up_rxns,'down_rxns':reg_analysis.down_rxns}
print(params_rxns)



# Code extracted from Jupyter Notebook cell

# here one can modify parameters of MiNEA 
path_to_params = '../input/task_enrichment_gem_params.yaml'



# Code extracted from Jupyter Notebook cell

# apply MiNEA with up and down regulated reactions
task_enrich = TaskEnrichment(cobra_model,path_to_params,params_rxns)

task_enrich.run()



# Code extracted from Jupyter Notebook cell

# run MiNEA for low and high expression reactions 
gene_exp={'gene_id':ery_df['Gene_id'].to_list(),'exp_val':ery_df['ERY_group_fpkm'].to_list(),'high_cutoff':0.15,'low_cutoff':0.15}

exp_analysis=ReactionExp(cobra_model,gene_exp=gene_exp)

params_rxns={'high_rxns':exp_analysis.high_rxns,'low_rxns':exp_analysis.low_rxns}

task_enrich = TaskEnrichment(cobra_model,path_to_params,params_rxns)

task_enrich.run()



# Code extracted from Jupyter Notebook cell



