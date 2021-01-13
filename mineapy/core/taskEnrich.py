from .enumTask import EnumerateTask
from .utils import remove_blocked_reactions, set_medium
import yaml
from .mhygepdf import mhygepdf
import statsmodels.stats.multitest as mtest
from scipy.stats import hypergeom
# from .logger import get_bistream_logger
from .numerics import BIGM,EPSILON

from cobra import Model
import pickle,os
import pandas as pd

from pytfa.optim.constraints import (
    SimultaneousUse,
    BackwardDirectionCoupling,
    ForwardDirectionCoupling,
)
from pytfa.optim.variables import (
    ForwardUseVariable,
    BackwardUseVariable,
)



from .thermo_model import ThermoModel_WithoutInfo

class TaskEnrichment():
    """
        A class to perform metabolic task enrichment analysis
    """

    def __init__(self,gem,parameters_path, params_rxn,inplace=False):

        # params_rxns={'up_rxns':reg_analysis.up_rxns,'down_rxns':reg_analysis.down_rxns, 'de_rxns':[de_rxns]}
        # params_rxns={'high_rxns':reg_analysis.up_rxns,'low_rxns':reg_analysis.down_rxns}
        self.read_parameters(parameters_path)
        # If inplace is True, no deepcopy is performed : the modifications are applied directly onto the gem
        prepared_gem = set_medium(gem, self.params['medium'], inplace)

        self._gem = prepared_gem
        self._params_rxn = params_rxn
        # This one is used to perform the lumping
        self._source_gem = prepared_gem.copy()
        ###


        ### test model id TFA or not
        if hasattr(self._source_gem , 'logger'):
            self.logger = self._gem.logger
            self._model_type='tfa'
        else:
            # self.logger = get_bistream_logger('FBA_model' + str(self._gem.id))
            ## we are going to add Forward and reverse use variable for the reactions
            tfa_model=ThermoModel_WithoutInfo(prepared_gem)

            tfa_model.convert_withoutInfo()
            self._source_gem = tfa_model
            self.logger = tfa_model.logger
            self._model_type='fba'

        ###
        self.fill_default_params()
        self.set_solver()


    def read_parameters(self, parameters_path):
        with open(parameters_path, 'r') as stream:
            try:
                self.params = yaml.safe_load(stream)
                print("Opened parameters file")
            except yaml.YAMLError as exc:
                print(exc)

    def fill_default_params(self):
        # If auto is activated, automatically extracts inorganics from the gem
        if "growth_rate" not in self.params or self.params["growth_rate"] == "auto":
            self.logger.info("Setting minimal growth rate to 95% of the TFA solution")
            obj_val = self._source_gem.slim_optimize()
            self.logger.info("Setting minimal growth rate to {}".format(obj_val))
            self.params["growth_rate"] = 0.95*obj_val
        if "force_solve" not in self.params:
            self.params["force_solve"] = False
        if "timeout" not in self.params:
            self.logger.info("Using default timeout : 3600s")
            self.params["timeout"] = 3600
        if "diverge_min" not in self.params:
            self.logger.info("Using default diverge_min : 1")
            self.params["diverge_min"] = 1
        if "feasibility" not in self.params:
            self.logger.info("Using default solver feasibility : 1e-9")
            self.params["feasibility"] = 1e-9
        else:
            # numbers like 1e-9 are detected as strings by yaml module
            # to enable their use, we cast them into floats
            try:
                self.params["feasibility"] = float(self.params["feasibility"])
            except ValueError as v:
                self.logger.error(v)

    def set_solver(self):
        if "solver" not in self.params or self.params["solver"].lower() == "auto":
            return None
        elif 'gurobi' in self.params["solver"].lower():
            solver = 'gurobi'
        elif 'cplex' in self.params["solver"].lower():
            solver = 'cplex'
        elif 'glpk' in self.params["solver"].lower():
            solver = 'glpk'
        else:
            solver = self.params["solver"]

        self._gem.solver = solver
        self._source_gem.solver = solver

    def get_enrichment_table_for_de(self,mins,de_rxns,de_col='deregulated'):
        ''' get enrichment of deregulated reactions (up or down regulated)'''
        mets_synthesis=[]
        alternatives=[]
        met_description=[]
        pvalues=[]
        probabilities=[]
        percent_de=[]
        rxns_num=[]
        list_de_in_min=[]
        list_rxn_in_min=[]
        gem=self._source_gem

        total_rxns=[]
        model_rxns=[]



        for rxn in gem.reactions:
            model_rxns.append(rxn.id)
            if rxn.gene_reaction_rule!="":
                total_rxns.append(rxn.id)
        total_rxns=set(total_rxns)
        model_rxns=set(model_rxns)
        without_gpr_rxns=[]

        for k,v in mins.items():
            #print(k,'\t',len(v[0]),len(v[1]),'\t',len(set(v[0].keys())-set(v[1].keys()))+len(set(v[1].keys())-set(v[0].keys())))

            for idx,a_min_dict in enumerate(v):
                mets_synthesis.append(k.id)
                # get metabolites
                try:
                    m=gem.metabolites.get_by_id(k.id)
                    met_description.append(m.name)
                except:
                    met_description.append(k)

                alternatives.append('alternative_%s'%str(idx+1))
                ## now we want to find up_regulated rxns and down_regulated reactions
                without_GPR=set([key.id for key in a_min_dict.keys()])
                rxns_in_min=set([key.id for key in a_min_dict.keys()]) & total_rxns
                de_in_min=rxns_in_min & set(de_rxns)
                percent_de.append(len(de_in_min)/len(rxns_in_min))
                rxns_num.append(len(rxns_in_min))
                list_de_in_min.append(','.join(de_in_min))
                list_rxn_in_min.append(','.join(rxns_in_min))

                without_gpr_rxns.append(len(without_GPR))
                ## total up and down reactions
                total_de_rxns=set(de_rxns)
                total_nochange_rxns=total_rxns-total_de_rxns
                ## apply hypergeometric function
                #M=len(total_rxns) ## number of gpr associated reactions (this can be true)
                M=len(model_rxns) ## number of gpr associated reactions
                n=len(de_rxns) ## total number of de regulated reactions
                N=len(rxns_in_min) ## number of reactions in a min
                x=len(de_in_min) ## number of deregulated in a min
                prob = hypergeom.cdf(x, M, n, N)

                probabilities.append(prob)
                ## calculate p-value
                pvalues.append(1-prob)


        df=pd.DataFrame()
        df['Metabolic task id']=mets_synthesis
        df['Alternative networks']=alternatives
        df['Metabolic task desciption']=met_description
        df[de_col+' rxns(%)']=percent_de
        df['GPR rxns']=rxns_num
        df['total # of rxns']=without_gpr_rxns
        df['p_values']=pvalues
        df['fdr_bh']=mtest.multipletests(pvalues, alpha=0.05, method='fdr_bh')[1]
        df[de_col+' rxns']=list_de_in_min
        df['rxns in a min']=list_rxn_in_min

        ## sort dataframe with p-values

        res_df=df.sort_values(by=['p_values'],ascending=True)
        res_df.to_csv(os.getcwd()+'/'+de_col+'_rxns.txt',sep='\t')



    def get_enrichment_table(self,mins,up_rxns,down_rxns,up_col='upregulated', down_col='downregulated'):
        ''' up and down can be high and low '''
        df=pd.DataFrame()
        mets_synthesis=[]
        alternatives=[]
        met_description=[]
        pvalues=[]
        probabilities=[]
        percent_up=[]
        percent_down=[]
        rxns_num=[]
        list_up_in_min=[]
        list_down_in_min=[]
        list_rxn_in_min=[]
        gem=self._source_gem
        total_rxns=[]
        model_rxns=[]
        for rxn in gem.reactions:
            model_rxns.append(rxn.id)
            if rxn.gene_reaction_rule!="":
                total_rxns.append(rxn.id)
        model_rxns=set(model_rxns)
        total_rxns=set(total_rxns)
        without_gpr_rxns=[]
        for k,v in mins.items():
            #print(k,'\t',len(v[0]),len(v[1]),'\t',len(set(v[0].keys())-set(v[1].keys()))+len(set(v[1].keys())-set(v[0].keys())))

            for idx,a_min_dict in enumerate(v):
                mets_synthesis.append(k.id)
                # get metabolites
                try:
                    m=gem.metabolites.get_by_id(k.id)
                    met_description.append(m.name)
                except:
                    met_description.append(k)

                alternatives.append('alternative_%s'%str(idx+1))
                ## now we want to find up_regulated rxns and down_regulated reactions
                without_GPR=set([key.id for key in a_min_dict.keys()])
                rxns_in_min=set([key.id for key in a_min_dict.keys()]) & total_rxns
                up_in_min=rxns_in_min & set(up_rxns)
                down_in_min=rxns_in_min & set(down_rxns)
                #nochange_in_min=rxns_in_min-up_in_min-down_in_min
                nochange_in_min=without_GPR-up_in_min-down_in_min
                # test division by error
                if len(rxns_in_min)>0:
                    percent_up.append((len(up_in_min)/len(rxns_in_min))*100)
                    percent_down.append((len(down_in_min)/len(rxns_in_min))*100)
                else:
                    percent_up.append(0)
                    percent_down.append(0)
                rxns_num.append(len(rxns_in_min))
                list_up_in_min.append(','.join(up_in_min))
                list_down_in_min.append(','.join(down_in_min))
                list_rxn_in_min.append(','.join(rxns_in_min))
                without_gpr_rxns.append(len(without_GPR))
                ## total up and down reactions
                total_up_rxns=set(up_rxns)
                total_down_rxns=set(down_rxns)
                #total_nochange_rxns=total_rxns-total_up_rxns-total_down_rxns ## this can be true
                total_nochange_rxns=total_rxns-total_up_rxns-total_down_rxns
                m=[len(total_up_rxns),len(total_down_rxns),len(total_nochange_rxns)];
                n=[len(up_in_min),len(down_in_min),len(nochange_in_min)];
                prob=mhygepdf(m,n);
                probabilities.append(prob)
                ## calculate p-value

                pval=0;
                for it1 in range(len(up_in_min),min(len(total_up_rxns),len(rxns_in_min))):
                    for it2 in range(len(down_in_min)):
                        if (it1+it2<len(rxns_in_min)):
                            pval=pval+mhygepdf(m,[it1,it2,len(rxns_in_min)-it1-it2]);
                pvalues.append(pval)
        df['Metabolic task id']=mets_synthesis
        df['Alternative networks']=alternatives
        df['Metabolic task desciption']=met_description
        df[up_col+' rxns(%)']=percent_up
        df[down_col+ ' rxns(%) ']=percent_down
        df['GPR rxns']=rxns_num
        df['total # of rxns']=without_gpr_rxns
        df['p_values']=pvalues
        df['fdr_bh']=mtest.multipletests(pvalues, alpha=0.05, method='fdr_bh')[1]
        df[up_col+' rxns']=list_up_in_min
        df[down_col+ ' rxns']=list_down_in_min
        df['rxns in a min']=list_rxn_in_min

        ## sort dataframe with p-values
        res_df=df.sort_values(by=['p_values'],ascending=True)
        res_df.to_csv(os.getcwd()+'/'+up_col+'_rxns.txt',sep='\t')





    def apply_enrichment(self,mins):
        ''' apply enrichment analysis '''
        params_rxn=self._params_rxn
        if ('up_rxns' in params_rxn.keys()) or ('down_rxns' in params_rxn.keys()):
            if (len(params_rxn['up_rxns'])>0) or (len(params_rxn['down_rxns'])>0):
                self.get_enrichment_table(mins,params_rxn['up_rxns'],params_rxn['down_rxns'],up_col='upregulated', down_col='downregulated')
                self.get_enrichment_table(mins,params_rxn['down_rxns'],params_rxn['up_rxns'],up_col='downregulated', down_col='upregulated')
                de_rxns=params_rxn['up_rxns']+params_rxn['down_rxns']
                self.get_enrichment_table_for_de(mins,de_rxns,de_col='deregulated')
            else:
                self.logger.info("empty up and down regulated reactions")
        elif ('high_rxns' in params_rxn.keys()) or ('low_rxns' in params_rxn.keys()):
            if (len(params_rxn['high_rxns'])>0) or (len(params_rxn['low_rxns'])>0):
                self.get_enrichment_table(mins,params_rxn['high_rxns'],params_rxn['low_rxns'],up_col='high_expressed', down_col='low_expressed')
            else:
                self.logger.info("empty up and down regulated reactions")
        else:
            self.logger.info("you need to reaction parameter for enrichment analysis 1) for condition change use up_ and down_rxns 2) for context use high_ and low_rxns")







    def run(self):
        # Extracting parameters

        biomass_rxn_ids = self.params["biomass_rxns"]
        biomass_rxns = [self._gem.reactions.get_by_id(x) for x in biomass_rxn_ids]
        main_bio_rxn = biomass_rxns[0]

        growth_rate = self.params["growth_rate"]
        cofactor_pairs = self.params["cofactor_pairs"]
        # Flatten cofactor_pairs list
        cofactors = [cofactor for pair in cofactor_pairs for cofactor in pair]
        lump_method = self.params["lump_method"]

        force_solve = self.params["force_solve"]
        timeout = self.params["timeout"]
        self._gem.solver.configuration.tolerances.feasibility = self.params["feasibility"]
        self._gem.solver.configuration.tolerances.integrality = self.params["feasibility"]
        self._source_gem.solver.configuration.tolerances.feasibility = self.params["feasibility"]
        self._source_gem.solver.configuration.tolerances.integrality = self.params["feasibility"]

        #### core reactions
        self.params['model_type']=self._model_type

        self.logger.info("Enumerating minmal networks ...")
        mt_model = EnumerateTask(self._source_gem, self.params) #### metabolic task model

        mins = mt_model.compute_mins(force_solve, method = lump_method)

        self.logger.info("Enumeration is done. now save the result in pickle file")
        pickle.dump(mins,open(os.getcwd()+'/mins.pickle','wb'))
        self.logger.info("Enumerating minimal networks ...")
        self.apply_enrichment(mins)
        # for k,v in mins.items():
        #     print(k,'\t',len(v[0]),len(v[1]),'\t',len(set(v[0].keys())-set(v[1].keys()))+len(set(v[1].keys())-set(v[0].keys())))
