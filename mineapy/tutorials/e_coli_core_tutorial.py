import mineapy
from cobra.io import load_matlab_model
import pandas as pd
from mineapy.core.taskEnrich import TaskEnrichment
from mineapy.core.thermo_model import ThermoModel_WithoutInfo
from mineapy.core.rxnExp import ReactionExp
cobra_model= load_matlab_model('./models/e_coli_core.mat')
genes=[g.id for g in cobra_model.genes]


path_to_params = './input/task_enrichment_params.yaml'

context_df=pd.read_csv('./input/context.txt',sep='\t')
condition_df=pd.read_csv('./input/condition.txt',sep='\t')
### working with condition comparisons
gene_reg={'gene_id':condition_df['geneid'].to_list(),'fold_change':condition_df['fold change'].to_list(),'up_cutoff':1.35,'down_cutoff':float(1/2.5)}

reg_analysis=ReactionExp(cobra_model,gene_reg=gene_reg)

gene_exp={'gene_id':context_df['geneid'].to_list(),'exp_val':context_df['exp_val'].to_list(),'high_cutoff':0.15,'low_cutoff':0.15}

exp_analysis=ReactionExp(cobra_model,gene_exp=gene_exp)

#params_rxns={'up_rxns':reg_analysis.up_rxns,'down_rxns':reg_analysis.down_rxns}
params_rxns={'high_rxns':exp_analysis.high_rxns,'low_rxns':exp_analysis.low_rxns}

task_enrich = TaskEnrichment(cobra_model,path_to_params,params_rxns)

task_enrich.run()
