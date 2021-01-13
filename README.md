## MiNEApy
This software can be used to perform minimal network enrichment analysis.

Vikash Pandey and Vassily Hatzimanikatis. Investigating the deregulation of metabolic tasks via Minimum Network Enrichment Analysis (MiNEA) as applied to nonalcoholic fatty liver disease using mouse and human omics data.DOI:https://doi.org/10.1371/journal.pcbi.1006760

## Requirements
You will need to have Git-LFS in order to properly download some binary files:

git clone https://github.com/EPFL-LCSB/pytfa.git /path/to/pytfa
git clone https://github.com/vpandey-om/mineapy.git /path/to/mineapy
cd /path/to/mineapy
python setup.py install

This module was developed in Python 3.7, and it is recommended to run Python 3.7 to run commercial solvers such as Gurobi and CPLEX.

This module requires COBRApy and pyTFA, as well as optlang to work properly. The installer should take care of that for you. You might also want to install a dedicated solver. GLPK, CPLEX and Gurobi are supported.

## Setup
This step is not required if you're using the container, which bundles all this.

You can install this module with pip:

For Python 3, you might have to use pip3 instead of pip

pip3 install MiNEApy
or from source

git clone https://github.com/vpandey-om/mineapy.git /path/to/mineapy
pip3 install -e /path/to/mineapy
Quick start
Two tutorial files detail thoroughly normal usages of the mineapy package. They can be found at:

mineapy
└── tutorials
    ├── e_coli_core.py
    ├── e_coli_gem.py
More information can be found here.



## Usage
First, create your COBRApy model. Make sure to define the additional values required by pyTFA, as said in the "Models" page of the documentation.

If you already have a Matlab model with thermodynamic data, you might want to use pytfa.io.import_matlab_model. Otherwise, have a look at the COBRApy documentation, then add the required properties.

If you're using a specific solver, don't forget to tell COBRApy about it by setting the solver property of your model to the name of your solver. See the COBRApy documentation for more information about this.

## Thermodynamic database
You also need a thermodynamic database. Use thermoDBconverter.py if you have a thermodynamic database from Matlab you wish to import to Python.

Thermodynamic databases are stored in .thermodb files and can be easily loaded with pytfa.io.load_thermoDB.

## Example script
Here is an example script :

## import required modules
import mineapy
from cobra.io import load_matlab_model
import pandas as pd
from mineapy.core.taskEnrich import TaskEnrichment
from mineapy.core.thermo_model import ThermoModel_WithoutInfo
from mineapy.core.rxnExp import ReactionExp
## load e_coli_core model
cobra_model= load_matlab_model('./models/e_coli_core.mat')
genes=[g.id for g in cobra_model.genes]

## Minea parameters
path_to_params = './input/task_enrichment_params.yaml'

## add condition- or context-specific expression data
context_df=pd.read_csv('./input/context.txt',sep='\t')
condition_df=pd.read_csv('./input/condition.txt',sep='\t')
## get genes that are regulated between different conditions
gene_reg={'gene_id':condition_df['geneid'].to_list(),'fold_change':condition_df['fold change'].to_list(),'up_cutoff':1.35,'down_cutoff':float(1/2.5)}

reg_analysis=ReactionExp(cobra_model,gene_reg=gene_reg)

## set cut off for example 15 % top and 15 % bottom in ranking  
gene_exp={'gene_id':context_df['geneid'].to_list(),'exp_val':context_df['exp_val'].to_list(),'high_cutoff':0.15,'low_cutoff':0.15}

exp_analysis=ReactionExp(cobra_model,gene_exp=gene_exp)


params_rxns={'high_rxns':exp_analysis.high_rxns,'low_rxns':exp_analysis.low_rxns}

## Apply enrichment algorithms
task_enrich = TaskEnrichment(cobra_model,path_to_params,params_rxns)

task_enrich.run()



The software in this repository is put under an APACHE-2.0 licensing scheme - please see the LICENSE file for more details.
