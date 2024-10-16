import mineapy
from cobra.io import load_matlab_model,load_json_model
import pandas as pd
from mineapy.core.taskEnrich import TaskEnrichment
from mineapy.core.thermo_model import ThermoModel_WithoutInfo
from mineapy.core.rxnExp import ReactionExp
from cobra.flux_analysis.variability import flux_variability_analysis
import sys
import os

# Add the path to the tutorials directory to sys.path
sys.path.append(os.path.abspath(os.path.join('..')))

#print(os.path.abspath(os.path.join('..')))
from data_utility import getErythromycin_data

## load the gene expression data for gene expression   Erythromycin 
ery_df=getErythromycin_data()

ery_df['significant(ERYvsETH)'].unique()
# Convert the 'log2FoldChange(ERYvsETH)' column to numeric, coercing errors to NaN
ery_df['log2FoldChange(ERYvsETH)'] = pd.to_numeric(ery_df['log2FoldChange(ERYvsETH)'], errors='coerce')

## find up and down regulated genes 
up_down_df=ery_df[(ery_df['significant(ERYvsETH)']=='UP')|(ery_df['significant(ERYvsETH)']=='DOWN')]
up_down_df['fold_change']=2 ** up_down_df['log2FoldChange(ERYvsETH)']

print(up_down_df)

# load core metabolic genome_scale_model
cobra_model = load_json_model('../models/iJO1366.json')


genes=[g.id for g in cobra_model.genes]


### working with condition comparisons
gene_reg={'gene_id':up_down_df['Gene_id'].to_list(),'fold_change':up_down_df['fold_change'].to_list(),'up_cutoff':1.1,'down_cutoff':0.9}

reg_analysis=ReactionExp(cobra_model,gene_reg=gene_reg)
reg_analysis
params_rxns={'up_rxns':reg_analysis.up_rxns,'down_rxns':reg_analysis.down_rxns}
print(params_rxns)

# here one can modify parameters of MiNEA 
path_to_params = '../input/task_enrichment_gem_params.yaml'

# apply MiNEA with up and down regulated reactions
task_enrich = TaskEnrichment(cobra_model,path_to_params,params_rxns)

task_enrich.run()