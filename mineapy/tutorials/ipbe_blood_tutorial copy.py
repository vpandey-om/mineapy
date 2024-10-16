# load e_coli core model
# test whetther they are FBA or TFA
# if not TFA then convert in TFA add Forward and reverse fluxes
# Perform flux variability analysis FBA and TFA
#
import sys
import pandas as pd

## add root folder of mineapy to sys.path folder
# for the moment it is hard coded
sys.path.append('/Users/vikash/python_Projects/mineapy/scripts')
print(sys.path)
##
import numpy as np
import pytfa
import pickle
from pytfa.analysis import  variability_analysis
from pytfa.io import import_matlab_model, load_thermoDB
from cobra.io import load_matlab_model
from core.taskEnrich import TaskEnrichment
from core.thermo_model import ThermoModel_WithoutInfo
from core.rxnExp import ReactionExp

data=pd.read_csv('/Users/vikash/python_Projects/mineapy/ipbe/male_female.txt',sep='\t')



ipbe_blood=pickle.load(open('/Users/vikash/python_Projects/mineapy/models/pbe_model/ipbeblood_py.pickle','rb'))
setattr(ipbe_blood, 'annotation',ipbe_blood._annotation)
for item in ipbe_blood.metabolites:
    setattr(item, 'annotation',item._annotation)
for item in ipbe_blood.reactions:
    setattr(item, 'annotation',item._annotation)
for item in ipbe_blood.genes:
    setattr(item, 'annotation',item._annotation)


#ipbe_blood2= pickle.load(open('/Users/vikash/Documents/Projects/MiNEApy/pytfa/models/small_ecoli_thermo.pickle','rb'))

# setattr(ipbe_blood, 'annotation', {})
#remove lower_bound of biomass to zero
biomass = ipbe_blood.reactions.get_by_id('biomass')
biomass.lower_bound=0
tfa_solution = ipbe_blood.optimize()
tfa_value = tfa_solution.objective_value
print('TFA Solution found : {0:.5g}'.format(tfa_value))

## get biomass reactions
biomass_rxn='biomass'

gene_ids=[g.id for g in ipbe_blood.genes]
values=[]

for g in gene_ids:
    if 'putative_' in g:
        g=g.replace('putative_','')
    tmp=data[data['GeneID']==g]
    if tmp.empty:
        values.append(np.nan)
    else:
        values.append(tmp['male_O'].values[0])

df=pd.DataFrame()
df['geneid']=gene_ids
df['exp_val']=values


path_to_params = '/Users/vikash/python_Projects/mineapy/ipbe/Minea_parameter_ipbe.yaml'


gene_exp={'gene_id':df['geneid'].to_list(),'exp_val':df['exp_val'].to_list(),'high_cutoff':0.15,'low_cutoff':0.15}

exp_analysis=ReactionExp(ipbe_blood,gene_exp=gene_exp)

#params_rxns={'up_rxns':reg_analysis.up_rxns,'down_rxns':reg_analysis.down_rxns}
params_rxns={'high_rxns':exp_analysis.high_rxns,'low_rxns':exp_analysis.low_rxns}

task_enrich = TaskEnrichment(ipbe_blood,path_to_params,params_rxns)

task_enrich.run()
