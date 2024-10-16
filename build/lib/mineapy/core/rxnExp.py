
from .logger import get_logger
import numpy as np
from scipy.stats import gmean
import pandas as pd
import math


class ReactionExp():
    """
        This class is used for reaction expression analysis
    """

    def __init__(self,gem,gene_reg={},gene_exp={}):
        ## gene_reg={'gene_id':[],'gene_fold_chnage':[],'up':value (default 2),'down':value (default 2)}
        ## gene_reg={'gene_id':[],'gene_expression':[],'top_percent':value (default 15%),'bottom_percent': value (default 15%)}
        analysis=self.test_type_of_analysis(gene_reg,gene_exp)
        if 'reg' in analysis.keys() and analysis['reg']=='yes':
            # apply reaction regulation
            up_rxns,down_rxns,rxn_df=self.getReactionRegulation(gem,gene_reg)
            self.up_rxns=up_rxns
            self.down_rxns=down_rxns
            self.rxn_df=rxn_df
        elif 'context' in analysis.keys() and analysis['context']=='yes':
            high_rxns,low_rxns,rxn_df=self.getReactionExpression(gem,gene_exp)
            self.high_rxns=high_rxns
            self.low_rxns=low_rxns
            self.rxn_df=rxn_df

    def test_type_of_analysis(self,gene_reg,gene_exp):
        ''' set parameters'''
        analysis={}
        logger = get_logger('expression_logger')
        self.logger=logger
        if bool(gene_reg):
            self.logger=logger.info('start analysis of comparison between conditions .....')
            ## now we want to check gene_id and values
            cond1=('gene_id' in gene_reg.keys())
            cond2=('fold_change' in gene_reg.keys())
            if  cond1 and  cond2 and len(gene_reg['gene_id'])>0 and len(gene_reg['fold_change'])>0:
                analysis={'reg':'yes'}
            else:
                self.logger.error('problem in input file please prepare as condition.txt')
        elif bool(gene_exp):
            self.logger.info('start analysis of context (tissue-specific).....')

            cond1=('gene_id' in gene_exp.keys())
            cond2=('exp_val' in gene_exp.keys())
            if  cond1 and  cond2 and len(gene_exp['gene_id'])>0 and len(gene_exp['exp_val'])>0:
                analysis={'context':'yes'}
            else:
                self.logger.error('problem in input file please prepare as context.txt')


        elif (bool(gene_reg)) and ( bool(gene_reg)):
            self.logger.info('Asking for both context (tissue-specific) and conditions' + \
            'comparison at the same time.....\n please pass only one (either gene_reg or gene_exp)')
        else:
            self.logger.info('You do not want to do anay analysis...\n( please pass only one (either gene_reg or gene_exp))')
        return analysis

    def getReactionRegulation(self,gem,gene_reg):
        ''' apply gene protein reaction rule for gene expression '''
        # get list of gpr and list of reactions ids
        # set by default parameters
        if 'up_cutoff' not in gene_reg.keys():
            gene_reg['up_cutoff']=2
        if 'down_cutoff' not in gene_reg.keys():
            gene_reg['down_cutoff']=2
        gpr_assoicated_reaction=[]
        gprs=[]

        for rxn in gem.reactions:
        	if rxn.gene_reaction_rule!="":
        		gpr_assoicated_reaction.append(rxn)
        		gprs.append(rxn.gene_reaction_rule)
        if len(gprs)>0:
            gene_reg_dict=dict(zip(gene_reg['gene_id'],gene_reg['fold_change']))
            reaction_regulation=self.evalGPR(gprs,gpr_assoicated_reaction,gene_reg_dict,OR_str='mean',AND_str='geomean')
            # get not nan values
            reaction_regulation={k.id: v for k, v in reaction_regulation.items() if pd.Series(v).notna().all()}
            reaction_df = pd.DataFrame(list(reaction_regulation.items()),columns = ['Reaction id','Fold change'])
            sorted_reaction_rgulation_df=reaction_df.sort_values(by=['Fold change'],ascending=False)
            if not (sorted_reaction_rgulation_df.empty):
                up_rxns,down_rxns=self.applyCutoffOnRxndf(sorted_reaction_rgulation_df,gene_reg['up_cutoff'],gene_reg['down_cutoff'],colname='Fold change')
                return up_rxns,down_rxns,sorted_reaction_rgulation_df
            else:
                self.logger.error('Can not map gpr to values. please check gpr rules ')

        else:
            self.logger.error('In the model there is no gene protein reaction rule is associated')


    def getReactionExpression(self,gem,gene_exp):
        ''' apply gene protein reaction rule for gene expression '''
        # get list of gpr and list of reactions ids
        # set by default parameters
        if 'high_cutoff' not in gene_exp.keys():
            gene_exp['high_cutoff']=.15
        if 'low_cutoff' not in gene_exp.keys():
            gene_exp['low_cutoff']=.15
        gpr_assoicated_reaction=[]
        gprs=[]

        for rxn in gem.reactions:
        	if rxn.gene_reaction_rule!="":
        		gpr_assoicated_reaction.append(rxn)
        		gprs.append(rxn.gene_reaction_rule)
        if len(gprs)>0:
            gene_exp_dict=dict(zip(gene_exp['gene_id'],gene_exp['exp_val']))
            reaction_regulation=self.evalGPR(gprs,gpr_assoicated_reaction,gene_exp_dict,OR_str='sum',AND_str='min')
            # get not nan values
            reaction_regulation={k.id: v for k, v in reaction_regulation.items() if pd.Series(v).notna().all()}
            reaction_df = pd.DataFrame(list(reaction_regulation.items()),columns = ['Reaction id','Expression value'])
            sorted_reaction_rgulation_df=reaction_df.sort_values(by=['Expression value'],ascending=False)
            if not (sorted_reaction_rgulation_df.empty):
                up_rxns,down_rxns=self.applyHighLowCutoff(sorted_reaction_rgulation_df,gene_exp['high_cutoff'],gene_exp['low_cutoff'],colname='Expression value')
                return up_rxns,down_rxns,sorted_reaction_rgulation_df
            else:
                self.logger.error('Can not map gpr to values. please check gpr rules ')

        else:
            self.logger.error('In the model there is no gene protein reaction rule is associated')


    def applyCutoffOnRxndf(self,rxn_df,up_cut,down_cut,colname):
        ''' apply cutoff on reaction regulations '''

        up_df=rxn_df[rxn_df[colname]>=up_cut]
        down_df=rxn_df[rxn_df[colname]<=down_cut]
        up_rxns=up_df['Reaction id'].to_list()
        down_rxns=down_df['Reaction id'].to_list()
        return up_rxns,down_rxns


    def applyHighLowCutoff(self,rxn_df,up_cut,down_cut,colname):
        ''' apply cutoff on reaction regulations '''
        up_idx=[rxn_df.index[i] for i in range(math.floor(rxn_df.shape[0]*up_cut))]

        down_idx=[rxn_df.index[rxn_df.shape[0]-1-i] for i in range(math.floor(rxn_df.shape[0]*down_cut))]

        up_df=rxn_df.loc[up_idx,:].copy()
        down_df=rxn_df.loc[down_idx,:].copy()
        up_rxns=up_df['Reaction id'].to_list()
        down_rxns=down_df['Reaction id'].to_list()
        return up_rxns,down_rxns


    def evalGPR(self,gpr_id,reaction_id,gene_reg_dict,OR_str='max',AND_str='min'):
        ''' This function is written for evaluate GPR. OR and AND is replaces by function such as add, min'''
        #### INPUT ###
        # gpr_id ## list of gpr rules
        # reaction_id ## list of reaction ids
        # gene_df is DataFrame: genes and columns as values
        # OR_str='max'(by default). you can choose 'mean','sum'
        # AND_str='min', (by default). you can choose 'geomean'
        #### OUTPUT ###

        # if ('genes' in gene_df.columns) and ('exp_vals' in gene_df.columns):
        #     gene_reg_dict=dict(zip(gene_df['genes'],gene_df['exp_vals']))
        # else:
        #     print(''' Please give input for genes and their values.
        #      "genes" and "exp_vals" should be header for the sample_gene_expression file (Which is tab delimited).
        #      The genes column contains model genes and the exp_vals column contains expression measurement ''')
        #     exit()

        i=0
        varying_reaction={}
        for ids in gpr_id:
            reaction=reaction_id[i]
            i=i+1
            ids=ids.strip("\n")
            or_part=ids.split(" or ")

            or_vals=[]
            for word in or_part:
                and_part=word.split(" and ")
                ### if length of and_part is more than 1 means (and) is the part of the string
                and_vals=[]
                for rv in and_part:
                    rv=rv.strip()
                    rv=rv.strip("(")
                    rv=rv.strip(")")
                    if rv in gene_reg_dict.keys() and (not (gene_reg_dict[rv]==np.nan)):
                        and_vals.append(gene_reg_dict[rv])

                if len(and_vals)>0:
                    if AND_str=='min':
                        or_vals.append(np.min(and_vals))
                    elif AND_str=='geomean':
                        or_vals.append(gmean(and_vals))
                    else:
                        print("Put 'min' or 'geomean' for  AND string ")
                else:
                    print(" gene value is not found or nan for the reaction = %s" %reaction.id)
                    continue

            if len(or_vals)>0:
                ## now we ar going to apply sum rules
                if OR_str=='max':
                    final_val=np.max(or_vals)
                elif OR_str=='mean':
                    final_val=np.mean(or_vals)
                elif OR_str=='sum':
                    final_val=np.sum(or_vals)
                else:
                    print("Put 'max' or 'mean' or 'sum' for  AND string ")

            else:
                print('why or_vals are empty')
                final_val=np.nan

            varying_reaction[reaction]=final_val
        return varying_reaction
