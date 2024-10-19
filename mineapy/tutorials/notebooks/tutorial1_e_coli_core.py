# import libraries 

import mineapy

from cobra.io import load_matlab_model

import pandas as pd

from mineapy.core.taskEnrich import TaskEnrichment

from mineapy.core.thermo_model import ThermoModel_WithoutInfo

from mineapy.core.rxnExp import ReactionExp

import sys

import os



# Add the path to the tutorials directory to sys.path

sys.path.append(os.path.abspath(os.path.join('..')))



#print(os.path.abspath(os.path.join('..')))

from data_utility import getErythromycin_data
# ## Data Analysis Using MiNEApy (Core Metabolic Model)

# 

# This tutorial demonstrates the application of Minimum Network Enrichment Analysis (MiNEA) using the MiNEApy package, utilizing transcriptomic data from a study on the response of Escherichia coli (E. coli) to erythromycin, an antibiotic ([Bie et al., 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10100721/)). The analysis is divided into two parts:

# 

# 1. [**Effect of Erythromycin on Gene Regulation:**](#effect-of-erythromycin-on-gene-regulation)

# Gene expression levels from erythromycin-treated cells are compared with control cells to identify upregulated and downregulated metabolic networks in the genome-scale model.

#    

# 2. [**FPKM Values and Network Enrichment:**](#fpkm-values-and-network-enrichment)  

# The second analysis uses FPKM (Fragments Per Kilobase of transcript per Million mapped reads) values from the erythromycin-treated group to identify networks enriched by highly expressed genes. This method minimizes the influence of low-expression genes, demonstrating how MiNEApy enriches networks based on transcriptional responses under antibiotic stress.

# 

# These analyses provide insights into gene expression and network-level responses to antibiotic stress, illustrating how MiNEApy can be applied to systems biology research.

# <a name="effect-of-erythromycin-on-gene-regulation"></a>

# ### 1. Effect of Erythromycin on Gene Regulation

# To achieve this analysis step by step:

# - Load the gene expression data to identify the up- and down-regulated genes.

# - Map the genes to reactions using Gene-Protein-Reaction (GPR) rules. Compute up- and down-regulated reactions based on the gene expression changes.

# - Compute minimal networks by tuning parameters in the MiNEApy tool. 

# - Combine the identified reactions and minimal networks to perform enrichment analysis and discover which pathways or networks are significantly impacted by the antibiotic treatment.
## load the gene expression data for gene expression   Erythromycin 

ery_df=getErythromycin_data()



ery_df['significant(ERYvsETH)'].unique()

# Convert the 'log2FoldChange(ERYvsETH)' column to numeric, coercing errors to NaN

ery_df['log2FoldChange(ERYvsETH)'] = pd.to_numeric(ery_df['log2FoldChange(ERYvsETH)'], errors='coerce')



## find up and down regulated genes 

up_down_df=ery_df[(ery_df['significant(ERYvsETH)']=='UP')|(ery_df['significant(ERYvsETH)']=='DOWN')]

up_down_df['fold_change']=2 ** up_down_df['log2FoldChange(ERYvsETH)']



print(up_down_df)

# load core metabolic model of e_coli

cobra_model= load_matlab_model('../models/e_coli_core.mat')

genes=[g.id for g in cobra_model.genes] # this is the gene list
### working with condition comparisons

gene_reg={'gene_id':up_down_df['Gene_id'].to_list(),'fold_change':up_down_df['fold_change'].to_list(),'up_cutoff':1.1,'down_cutoff':0.9}



reg_analysis=ReactionExp(cobra_model,gene_reg=gene_reg)

reg_analysis

params_rxns={'up_rxns':reg_analysis.up_rxns,'down_rxns':reg_analysis.down_rxns}

print(params_rxns)
# ### Configuring the `parameter.yaml` for MiNEA

# 

# 

# In this section, the process of configuring the `parameter.yaml` file for running MiNEA (Minimum Network Enrichment Analysis) is described. This file controls the behavior of the analysis and enables adjustments based on specific research requirements. Below is an example `parameter.yaml` configuration, with explanations for each parameters:

# 

# 

# ```yaml

# # Biomass reactions: Specify the biomass reaction(s) used in your model.

# # Uncomment the lines below to define specific biomass reactions.

# biomass_rxns:

#     #- Biomass_Ecoli_core

#     #- Ec_biomass_iJO1366_WT_53p95M

#     - BIOMASS_Ecoli_core_w_GAM

# 

# # Define the minimal growth rate required. Setting to 'auto' uses 95% 

# # of the maximum growth rate determined by TFA or FBA.

# growth_rate: auto

# 

# # Define the Metabolic Tasks (MTs). These are specific metabolites 

# # for which you want to find the minimal networks.

# MTs: ['pyr_c', 'oaa_c','succ_c']

# 

# # Use Building Blocks (BBB) as Metabolic Tasks if set to 'yes'. This computes minimal networks for each substrate in the biomass reaction.

# MT_BBB: yes

# 

# # Define the medium composition. Specify which external metabolites 

# # are available in the medium. Reactions not listed will retain 

# # their model-defined bounds.

# medium:

#     #EX_glc__D_e: -10

# 

# # Cofactor pairs: Specify pairs of cofactors that are to be balanced 

# # during the analysis.

# cofactor_pairs:

#     - ['atp_c', 'adp_c']

#     - ['nad_c', 'nadh_c']

#     - ['nadp_c', 'nadph_c']

#     - ['accoa_c', 'coa_c']

# 

# # Lump method: This parameter will alllow maximum number of network size. 

# # Options include 'Min', 'Min+1', 'Min+30'.

# lump_method: Min+30

# 

# # Maximum number of lumped reactions per metabolite or metabolic task. Only applicable 

# # if 'Min' or 'Min+p' is selected as lump method. 

# # It defines the number of alternative networks that can be obtained for each metabolite.

# max_lumps_per_BBB: 2

# 

# # Diverge minimum: Specifies how much the solutions are allowed to deviate from the base solution. 

# # A higher value permits more variation between solutions, indicating the number of reactions that can differ between consecutive solutions. publication: Rezola A, Pey J, de Figueiredo LF, Podhorski A, Schuster S, Rubio A, Planes FJ. Selection of human tissue-specific elementary flux modes using gene expression data. Bioinformatics. 2013 Aug 15;29(16):2009-16. doi: 10.1093/bioinformatics/btt328. Epub 2013 Jun 6. PMID: 23742984.

# diverge_min: 5  # default=1 

# 

# # Solver options: Define which solver to use. 'auto' lets the program decide. 'gurobi', 'cplex','glpk' are the options. 

# solver: auto

# 

# # Force solving even if some tasks are infeasible.

# force_solve: False

# 

# # Timeout: Maximum time (in minutes) allocated for solving each problem.

# timeout: 300

# 

# # Feasibility tolerance: Define the precision for feasibility checks.

# feasibility: 1e-9

# 

# # Constraint method: Choose how to handle constraints. Options 

# # are 'flux', 'integer', or 'both'.

# constraint_method: both

# Modify parameters of MiNEA 

path_to_params = '../input/task_enrichment_params.yaml'
# apply MiNEA with up and down regulated reactions

task_enrich = TaskEnrichment(cobra_model,path_to_params,params_rxns)



task_enrich.run()

# ### Output Files Generated by MiNEA with up and down regulated reactions

# 

# Upon completing the Minimum Network Enrichment Analysis (MiNEA),  

# the algorithm will generate five key output files in the local directory  

# where the Jupyter notebook is executed. These files provide critical insights  

# into the metabolic tasks under study:

# 

# 1. **fva.pickle**:  

#    This file stores the results of the **Flux Variability Analysis (FVA)**  

#    conducted on the metabolic model. FVA helps determine the range of possible  

#    fluxes through each reaction under the given constraints, providing insight  

#    into the model's flexibility and possible reaction behaviors.

# 

# 2. **mins.pickle**:  

#    This file contains the **enumerated minimal networks** for each task.  

#    These minimal networks represent the smallest subsets of reactions that  

#    can carry out a specific metabolic task, offering a simplified view of  

#    the complex metabolic processes.

# 

# 3. **upregulated_rxns.txt**:  

#    This file contains the enrichment result tables for upregulated reactions.  

#    Each alternative minimal network is ranked based on p-values. The computation  

#    of p-values considers as many upregulated reactions as possible while  

#    minimizing the inclusion of downregulated reactions.

# 

# 4. **downregulated_rxns.txt**:  

#    This file contains the enrichment result tables for downregulated reactions.  

#    Similar to the upregulated reactions file, each alternative minimal network  

#    is ranked based on p-values. Here, the computation of p-values considers  

#    as many downregulated reactions as possible while minimizing the inclusion  

#    of upregulated reactions.

# 

# 5. **deregulated_rxns.txt**:  

#    This file contains the enrichment result tables for deregulated reactions.  

#    In this analysis, both upregulated and downregulated reactions are considered  

#    together without focusing on the direction of regulation. The p-values reflect  

#    the inclusion of both reaction types, offering a comprehensive view of  

#    the deregulation occurring in the metabolic network.

# <a name="fpkm-values-and-network-enrichment"></a>

# ### 2. FPKM Values and Network Enrichment

# To achieve this analysis step by step:

# - Load the gene expression data to identify the high- and low-expressed gene genes.

# - Using Gene-Protein-Reaction (GPR) rules, compute high- and low-expressed reactions based on the gene expression.

# - Compute minimal networks by tuning parameters in the MiNEApy tool. 

# - Combine high- and low-expressed reactions  and minimal networks to perform enrichment analysis and discover which pathways or networks are significantly enriched by high-expressed reactions.

# - publication:

# 1. Rezola A, Pey J, de Figueiredo LF, Podhorski A, Schuster S, Rubio A, Planes FJ. Selection of human tissue-specific elementary flux modes using gene expression data. Bioinformatics. 2013 Aug 15;29(16):2009-16. doi: 10.1093/bioinformatics/btt328. Epub 2013 Jun 6. PMID: 23742984.

# 2. Pandey V, Hatzimanikatis V. Investigating the deregulation of metabolic tasks via Minimum Network Enrichment Analysis (MiNEA) as applied to nonalcoholic fatty liver disease using mouse and human omics data. PLoS Comput Biol. 2019 Apr 19;15(4):e1006760. doi: 10.1371/journal.pcbi.1006760. PMID: 31002661; PMCID: PMC6493771. 
# run MiNEA for low and high expression reactions 

# A value of high_cutoff=0.15 indicates that the top 15% of highly expressed genes.

# and the bottom 15% of lowly expressed genes will be selected when low_cutoff=0.15.

gene_exp={'gene_id':ery_df['Gene_id'].to_list(),'exp_val':ery_df['ERY_group_fpkm'].to_list(),'high_cutoff':0.15,'low_cutoff':0.15}



exp_analysis=ReactionExp(cobra_model,gene_exp=gene_exp)



params_rxns={'high_rxns':exp_analysis.high_rxns,'low_rxns':exp_analysis.low_rxns}



task_enrich = TaskEnrichment(cobra_model,path_to_params,params_rxns)



task_enrich.run()
# ### Output Files Generated by MiNEA with High and Low Expressed Reactions

# 

# In addition to the output files previously described 

# (fva.pickle and mins.pickle), which store the results of 

# Flux Variability Analysis (FVA) and enumerated minimal networks 

# respectively, MiNEA generates a third file focused on 

# gene expression levels:

# 

# 3. **high_expressed_rxns.txt**:

# 

#    This file contains the enrichment result tables 

#    for highly expressed reactions. Each alternative 

#    minimal network is ranked based on p-values. The 

#    computation of p-values in this context prioritizes 

#    the inclusion of as many high expressed reactions 

#    as possible while minimizing the inclusion of low 

#    expressed reactions. This type of network analysis 

#    is particularly useful when identifying tissue-specific 

#    metabolic networks or when trying to understand the 

#    metabolic context of specific biological conditions.

