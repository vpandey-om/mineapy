from pytfa.analysis import  variability_analysis
from cobra import Reaction
#from cobra.flux_analysis.variability import flux_variability_analysis
from pytfa.optim.utils import symbol_sum
from pytfa.thermo.utils import is_exchange, check_transport_reaction
# from .utils import trim_epsilon_mets

from pytfa.optim.variables import ReactionVariable, BinaryVariable, get_binary_type
from pytfa.optim.constraints import ReactionConstraint, ForbiddenProfile

from numpy import sum, round,Inf

from optlang.interface import INFEASIBLE, TIME_LIMIT, OPTIMAL

from tqdm import tqdm

import pickle,os
import numpy as np

from collections import defaultdict, namedtuple
from .thermo_model import ThermoModel_WithoutInfo
from ..optim.constraints import SumFlux
from ..optim.variables import FluxSumVar

CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'


# Transforms (OnePerBBB --> oneperbbb), (one_per_bbb --> oneperbbb), etc ...
disambiguate = lambda s:s.lower().replace('_','')

Lump = namedtuple('Lump', ['id_', 'metabolites', 'subnetwork', 'gene_reaction_rule'])

class InfeasibleExcept(Exception):
    def __init__(self, status, feasibility):
        self.status = status
        self.feasibility = feasibility


class TimeoutExcept(Exception):
    def __init__(self, time_limit):
        self.time_limit = time_limit


class FluxKO(ReactionVariable, BinaryVariable):
    prefix = 'KO_'

    def __init__(self, reaction, **kwargs):
        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

# Define a new constraint type:
class UseOrKOInt(ReactionConstraint):
    prefix = 'UKI_'
# Define a new constraint type:
class UseOrKOFlux(ReactionConstraint):
    prefix = 'UKF_'



class EnumerateTask:
    """
        A class encapsulating the MiNEA algorithm
    """
    def __init__(self, tfa_model, params):
        """
        :type tfa_model: pytfa model

        :param biomass_rxns: list of biomass reactions
        :type biomass_rxns: [GEM.biomass_rxn.id]

        :param core_subsystems: list of Metabolic Tasks or BBB
        :type core_subsystems: [string]

        :param growth_rate: theoretical maximum specific growth rate in 1/hr units
        :type growth_rate: float

        :param timeout_limit: the maximum amount of time allowed to compute each optimization. Default is 3600s (1 hour)
        :type timeout_limit: float (seconds)
        """

        self._tfa_model = tfa_model

        self._param_dict = params
        self.init_params()
        # Set containing every Metabolic tasks
        # self._MTs = self.mts
        # Set containing every BBB reaction
        self._rBBB = list()
        # For each reaction

        if self.mt_bbb:
            for rxn in self._tfa_model.reactions:
                # If it's a BBB reaction
                #print(rxn.id,rxn.name)
                if rxn.id in self.biomass_rxns:
                    self._rBBB.append(rxn)
        # Growth rate
        self._growth_rate = self.growth_rate

        # TODO : solver choice
        # TODO default : solver du modele
        #self._solver = 'optlang-cplex'

        self._tfa_model.solver.configuration.timeout = self.timeout_limit
        print("Timeout limit is {}s".format(self.timeout_limit))

        try:
            fva=pickle.load(open(os.getcwd()+'/fva.pickle','rb'))
        except:
            self._tfa_model.logger.info(" Removing blocked reaction (that will help in optimization) ...")
            fva=self._get_blocked_reaction() # these are the blocked reactions

            pickle.dump(fva,open(os.getcwd()+'/fva.pickle','wb'))
        epsilon = self._tfa_model.solver.configuration.tolerances.feasibility
        blocked_df = fva[ (fva.max(axis=1).abs()<1*epsilon)
                & (fva.min(axis=1).abs()<1*epsilon)]
        self._blocked_reactions=blocked_df
        self._fva=fva

        # non_blocked_reactions
        not_blocked_rxns=[]
        for rxn in self._tfa_model.reactions:
            if rxn.id not in blocked_df.index.to_list():
                not_blocked_rxns.append(rxn)
        self._not_blocked_rxns=not_blocked_rxns

        ### remove bloked reactions
        ### initialize lb and ub from thermodynamic variability
        ### make FBA model

        # self._prune_and_make_FBA_model()
        # self._sinks=self._prepare_Mts()

        ## apply sum of fluxes
        # self._apply_sum_of_flux()
        #
        #
        #
        #
        #  # lumpgem binary variables to deactivate non-core reactions. The reaction is deactivated when the value of
        # # the variable is 1
        self._activation_vars = {rxn: self._tfa_model.add_variable(kind=FluxKO,
                                                                   hook=rxn,
                                                                   lb=0,
                                                                   ub=1,
                                                                   queue=False)
                                 for rxn in not_blocked_rxns}



        self._generate_usage_constraints()
        self._generate_objective()
        self._sinks=self._prepare_Mts()



    def _prune_FBA_model(self,tmp_model,bloked_rxn_ids):
        """
        Remove bloked reactions and make FBA model.
        working on reduced FBA model
        """
        tmp_model.remove_reactions(bloked_rxn_ids)
        return tmp_model


    def _prune_and_make_FBA_model(self):
        """
        Remove bloked reactions and make FBA model.
        working on reduced FBA model
        """
        #epsilon = self._tfa_model.solver.configuration.tolerances.feasibility
        rid_to_rm=self._blocked_reactions.index ## remove these reactions
        tmp_model=self._tfa_model.copy()
        tmp_model.remove_reactions(rid_to_rm)
        rid=[]
        for rxn in tmp_model.reactions:
            sink = tmp_model.reactions.get_by_id(rxn.id)
            # Activate reaction by setting its lower bound
            try:
                lb=self._fva.loc[rxn.id,'minimum']
                ub=self._fva.loc[rxn.id,'maximum']
            # get in fva
                sink.lower_bound=lb
                sink.upper_bound=ub
                rid.append(rxn.id)
            except:
                continue
        ## build FBA model


        tfa_model=ThermoModel_WithoutInfo(tmp_model)

        tfa_model.convert_withoutInfo()
        self._source_gem = tfa_model
        self._tfa_model = tfa_model
        self.logger = tfa_model.logger
        self.model_type='fba'
        self.orig_rid=rid


    def _apply_sum_of_flux(self):
        """
        we apply sum of flux as a variable and do flux variability analysis using optimum value as constraint
        """
        #epsilon = self._tfa_model.solver.configuration.tolerances.feasibility
        # get forward use and backward use variables
        forward_variables = [rxn.forward_variable  for rxn in self._source_gem.reactions]
        reverse_variables = [ rxn.reverse_variable for rxn in self._source_gem.reactions]
        #self._source_gem.reactions[0].forward_variable + self._source_gem.reactions[0].reverse_variable

        self._source_gem.add_variable(kind=FluxSumVar,hook=self._source_gem,id_='FluxsumVar',lb=0,ub=Inf)
        self._FluxsumVar=self._source_gem.variables[-1]
        expr = symbol_sum(forward_variables+reverse_variables)
        new_expr=expr-self._FluxsumVar
        self._source_gem.add_constraint(kind=SumFlux,hook = self._source_gem,id_ = 'SumFlux',expr = new_expr,lb=0,ub=0)

        self._source_gem.objective = self._source_gem.problem.Objective(self._FluxsumVar, direction='min')
        solution=self._source_gem.optimize()
        self._solution=solution
        if solution.status=='optimal':
            self._FluxsumVar.lb=solution.objective_value
            self._FluxsumVar.ub=solution.objective_value+0.1*solution.objective_value

    def _solve_sum_of_flux(self):
        """
        we apply sum of flux as a variable and do flux variability analysis using optimum value as constraint
        """
        #epsilon = self._tfa_model.solver.configuration.tolerances.feasibility
        # get forward use and backward use variables
        forward_variables = [rxn.forward_variable  for rxn in self._tfa_model.reactions]
        reverse_variables = [ rxn.reverse_variable for rxn in self._tfa_model.reactions]
        #self._source_gem.reactions[0].forward_variable + self._source_gem.reactions[0].reverse_variable

        self._tfa_model.add_variable(kind=FluxSumVar,hook=self._tfa_model,id_='FluxsumVar',lb=0,ub=Inf)
        self._FluxsumVar=self._tfa_model.variables[-1]
        expr = symbol_sum(forward_variables+reverse_variables)
        new_expr=expr-self._FluxsumVar
        self._tfa_model.add_constraint(kind=SumFlux,hook = self._tfa_model,id_ = 'SumFlux',expr = new_expr,lb=0,ub=0)

        self._tfa_model.objective = self._tfa_model.problem.Objective(self._FluxsumVar, direction='min')



    def _generate_usage_constraints(self):
        """
        Generate carbon intake related constraints for each non-core reaction
        For each reaction rxn : rxn.forward_variable + rxn.reverse_variable + activation_var * C_uptake < C_uptake
        """
        flux_methods = ['flux', 'fluxes', 'both']
        int_methods = ['int', 'integer', 'both']

        if self.constraint_method.lower() not in flux_methods + int_methods:
            raise ArgumentError('{} is not a correct constraint method. '
                                'Choose among [Flux, Integer, Both]. '
                                'If you do not know what to choose, go for Flux.'
                                'If it is too slow, go for integer.'
                                'If you get strange lumps, go for both'
                                .format(self.constraint_method))

        for rxn in self._not_blocked_rxns:
            activation_var = self._activation_vars[rxn]
            if self.constraint_method.lower() in flux_methods:
                bigM = 100
                reac_var = rxn.forward_variable + rxn.reverse_variable + activation_var * bigM
                # adding the constraint to the model
                self._tfa_model.add_constraint(kind=UseOrKOFlux,
                                               hook=rxn,
                                               expr=reac_var,
                                               ub=bigM,
                                               lb=0,
                                               queue=True)
            if self.constraint_method.lower() in int_methods:
                fu = self._tfa_model.forward_use_variable .get_by_id(rxn.id)
                bu = self._tfa_model.backward_use_variable.get_by_id(rxn.id)
                reac_var = fu + bu + activation_var
                # adding the constraint to the model
                self._tfa_model.add_constraint(kind=UseOrKOInt,
                                               hook=rxn,
                                               expr=reac_var,
                                               ub=1,
                                               lb=0,
                                               queue=True)


        # push constraints in one bulk (faster)
        self._tfa_model._push_queue()
        # refresh constraint fields
        self._tfa_model.repair()



    def _get_blocked_reaction(self):
        """
        We compute flux variability analysis. so that we could block those reactions
        """


        fva = variability_analysis(self._tfa_model, kind='reactions')



        return fva

    def init_params(self):
        self.biomass_rxns = self._param_dict["biomass_rxns"]

        self.growth_rate = self._param_dict["growth_rate"]

        self.timeout_limit = self._param_dict["timeout"]
        self.mts = self._param_dict["MTs"]
        self.mt_bbb = self._param_dict["MT_BBB"]
        self.constraint_method = self._param_dict["constraint_method"]
        self.model_type = self._param_dict["model_type"]
        self.cofactor_pairs = self._param_dict["cofactor_pairs"]
        self.diverge_min = self._param_dict["diverge_min"]

    def get_cofactor_adjusted_stoich(self,rxn):
        stoich_dict = {x.id:v for x,v in rxn.metabolites.items()}

        for a,b in self.cofactor_pairs:
            try:
                na = stoich_dict[a] # looks like -54 atp_c
                nb = stoich_dict[b] # looks like +53 adp_c

                n = na+nb # looks like -1

                if n == 0:
                    self._tfa_model.logger.warn(
                        'Cofactor pair {}/{} is equimolar in reaction {}'
                        .format(a,b,rxn.id))
                elif n > 0:
                    n = -n
                    self._tfa_model.logger.warn(
                        'Cofactor pair {}/{} looks inverted in reaction {}'
                        .format(a,b,rxn.id))

                stoich_dict[a] =  n # looks like 1
                stoich_dict[b] = -n # looks like -1
            except KeyError:
                pass
        return stoich_dict



    def _prepare_Mts(self):
        """
        For each BBB (reactant of the biomass reactions), generate a sink, i.e an unbalanced reaction BBB ->
        of which purpose is to enable the BBB to be output of the GEM
        :return: the dict {BBB: sink} containing every BBB (keys) and their associated sinks
        """
        all_sinks = {}
        print("Preparing metabolic tasks...")

        for bio_rxn in self._rBBB:
            stoich_dict = self.get_cofactor_adjusted_stoich(bio_rxn)
            for met in bio_rxn.reactants:
                stoech_coeff = stoich_dict[met.id]
                # stoech_coeff < 0 indicates that the metabolite is a reactant
                if (stoech_coeff < 0) and (met not in all_sinks.keys()):
                    sink = Reaction("MT_" + bio_rxn.id + "_" + met.id)
                    sink.name = "MT_" + bio_rxn.name + "_" + met.name
                    # Subsystem specific to BBB sinks
                    sink.subsystem = "MTs"

                    # A sink is simply a reaction which consumes the BBB
                    sink.add_metabolites({met: -1})

                    # The sinks will be activated later (cf compute_lumps), one at a time
                    # sink.knock_out()

                    # The stoechiometric coefficients will be used to define the lower bound of the sink,
                    # thus it must be stored
                    all_sinks[met] = (sink.id, -stoech_coeff)
                    self._tfa_model.add_reactions([sink])

                # reactant already seen
                elif stoech_coeff < 0:
                    # The BBB has already been associated to a sink, so we simply increase the bound of the sink
                    all_sinks[met][1] -= stoech_coeff

        # Must be called before changing the reaction.thermo['computed'] values

        for met in self._tfa_model.metabolites:
            ### get metabolite by id

            if met.id in self.mts:
                sink = Reaction("MT_"+ met.id)
                sink.name = "MT_" + met.name
                # Subsystem specific to BBB sinks
                sink.subsystem = "MTs"

                # A sink is simply a reaction which consumes the BBB
                sink.add_metabolites({met: -1})
                all_sinks[met] = (sink.id, 1)
                self._tfa_model.add_reactions([sink])

        if self.model_type=='tfa':
            self._tfa_model.prepare()
            # Must be called before optimization
            self._tfa_model.convert()
        else:
            self._tfa_model.convert_withoutInfo()


        # self._activation_vars = {rxn: self._tfa_model.add_variable(kind=FluxKO,
        #                                                            hook=rxn,
        #                                                            lb=0,
        #                                                            ub=1,
        #                                                            queue=False)
        #                          for rxn in self._tfa_model.reactions}
        # for ncrxn in self._rncore:
        #     ncrxn.thermo['computed'] = False

        return all_sinks


    def _generate_objective(self):
        """
        Generate and add the maximization objective : set as many activation variables as possible to 1
        When an activation variable is set to 1, the corresponding non-core reaction is deactivated
        """
        # Sum of binary variables to be maximized

        objective_sum = symbol_sum(list(self._activation_vars.values()))
        # Set the sum as the objective function
        self._tfa_model.objective = self._tfa_model.problem.Objective(objective_sum, direction='max')

    def compute_mins(self, force_solve=False, method='OnePerBBB'):
        """
        For each BBB (reactant of the biomass reaction), add the corresponding sink to the model, then optimize and
        lump the result into one lumped reaction
        :param force_solve: Indicates whether the computations must continue when one lumping yields a status "infeasible"
        :return: The dict {BBB: lump} containing every lumped reactions, associated to their BBBs
        """

        # # Must be called before optimization
        # self._tfa_model.convert()
        # self._tfa_model.objective_direction = 'min'

        epsilon = self._tfa_model.solver.configuration.tolerances.feasibility

        the_method = disambiguate(method)
        print('Min network method detected: {}'.format(the_method))

        # dict: {metabolite: lumped_reaction}
        lumps = {}

        self._tfa_model.objective_direction = 'max'

        sink_iter = tqdm(self._sinks.items(), desc = 'met')
        # import pdb; pdb.set_trace()
        for met_BBB, (sink_id, stoech_coeff) in sink_iter:

            # Cute stuff

            sink_iter.set_description('met={}'.format(met_BBB.id[:10]))
            sink_iter.refresh()

            sink = self._tfa_model.reactions.get_by_id(sink_id)
            # Activate reaction by setting its lower bound
            prev_lb = sink.lower_bound
            # min_prod = self._growth_rate * stoech_coeff
            # sink.lower_bound = min_prod - epsilon
            ### we are going to optimize
            self._tfa_model.objective = self._tfa_model.problem.Objective(sink.forward_variable, direction='max')
            n_da = self._tfa_model.slim_optimize()

            try:
                # Timeout reached
                if self._tfa_model.solver.status == TIME_LIMIT:
                    raise TimeoutExcept(self._tfa_model.solver.configuration.timeout)
                # Not optimal status -> infeasible
                elif self._tfa_model.solver.status != OPTIMAL:
                    raise InfeasibleExcept(self._tfa_model.solver.status,
                                           self._tfa_model.solver.configuration.tolerances.feasibility)

            except (TimeoutExcept, InfeasibleExcept) as err:
                print('Can not Produced {}'.format(sink_id))
                self.lumps = lumps
                return lumps

            sink.lower_bound = sink.flux*0.8
            self._generate_objective()
            if the_method == 'oneperbbb':
                this_lump = self._lump_one_per_bbb(met_BBB, sink, force_solve)
                lumped_reactions = [this_lump] if this_lump is not None else list()
            elif the_method.startswith('min+'):
                try:
                    p = int(the_method.replace('min+',''))
                except ValueError:
                    raise ValueError('Min+p method must have p as an integer')
                lumped_reactions = self._lump_min_plus_p(met_BBB, sink, p, force_solve)
            elif the_method.startswith('min'):
                lumped_reactions = self._lump_min_plus_p(met_BBB, sink, 0, force_solve)
            else:
                raise ValueError('Lumping method not recognized: {}. '
                                 'Valid methods are '
                                 'OnePerBBB, Min, Min+p, p natural integer'
                                 .format(the_method))


            if not lumped_reactions:
                continue

            lumps[met_BBB] = lumped_reactions
            # Deactivating reaction by setting both bounds to 0
            sink.lower_bound = prev_lb
            # sink.knock_out()

        self.lumps = lumps
        return lumps

    def _lump_one_per_bbb_mix(self, met_BBB, sink, force_solve):
        """

        :param met_BBB:
        :param sink:
        :param force_solve:
        :return:
        """
        n_da = self._tfa_model.slim_optimize()

        try:
            # Timeout reached
            if self._tfa_model.solver.status == TIME_LIMIT:
                raise TimeoutExcept(self._tfa_model.solver.configuration.timeout)
            # Not optimal status -> infeasible
            elif self._tfa_model.solver.status != OPTIMAL:
                raise InfeasibleExcept(self._tfa_model.solver.status,
                                       self._tfa_model.solver.configuration.tolerances.feasibility)

        except (TimeoutExcept, InfeasibleExcept) as err:
            # If the user want to continue anyway, suits him
            if force_solve:
                # Raise a warning

                self._solve_sum_of_flux()
                n_da = self._tfa_model.slim_optimize()
                try:
                    # Timeout reached
                    if self._tfa_model.solver.status == TIME_LIMIT:
                        raise TimeoutExcept(self._tfa_model.solver.configuration.timeout)
                        # Not optimal status -> infeasible
                    elif self._tfa_model.solver.status != OPTIMAL:
                        raise InfeasibleExcept(self._tfa_model.solver.status,
                                           self._tfa_model.solver.configuration.tolerances.feasibility)
                except (TimeoutExcept, InfeasibleExcept) as err:
                    return None

                lumped_reaction = self._build_lump_sum_flux(met_BBB, sink)
                return [lumped_reaction,'sum_flux']


            else:
                raise err

        print('Produced {}'.format(sink.flux),
              'with {0:.0f} reactions deactivated'.format(n_da))
        lumped_reaction = self._build_lump(met_BBB, sink)

        return lumped_reaction



    def _lump_one_per_bbb(self, met_BBB, sink, force_solve):
        """

        :param met_BBB:
        :param sink:
        :param force_solve:
        :return:
        """

        #n_da = self._tfa_model.slim_optimize()
        solution = self._tfa_model.optimize()
        n_da=solution.objective_value

        try:
            # Timeout reached
            if self._tfa_model.solver.status == TIME_LIMIT:
                raise TimeoutExcept(self._tfa_model.solver.configuration.timeout)
            # Not optimal status -> infeasible
            elif self._tfa_model.solver.status != OPTIMAL:
                raise InfeasibleExcept(self._tfa_model.solver.status,
                                       self._tfa_model.solver.configuration.tolerances.feasibility)

        except (InfeasibleExcept) as err:
            # If the user want to continue anyway, suits him
            if force_solve:
                # Raise a warning
                self._solve_sum_of_flux()
                n_da = self._tfa_model.slim_optimize()
                try:
                    # Timeout reached
                    if self._tfa_model.solver.status == TIME_LIMIT:
                        raise TimeoutExcept(self._tfa_model.solver.configuration.timeout)
                        # Not optimal status -> infeasible
                    elif self._tfa_model.solver.status != OPTIMAL:
                        raise InfeasibleExcept(self._tfa_model.solver.status,
                                           self._tfa_model.solver.configuration.tolerances.feasibility)
                except (TimeoutExcept, InfeasibleExcept) as err:
                    return None

                lumped_reaction = self._build_lump_sum_flux(met_BBB, sink)
                return [lumped_reaction,'sum_flux']


            else:
                raise err

        print('Produced {}'.format(sink.flux),
              'with {0:.0f} reactions deactivated'.format(n_da))
        lumped_reaction = self._build_lump(met_BBB, sink)

        return lumped_reaction

    def _build_lump(self, met_BBB, sink):
        """
        This function uses the current solution of self._tfa_model

        :param met_BBB:
        :param sink:
        :return:
        """

        epsilon_int = self._tfa_model.solver.configuration.tolerances.integrality
        epsilon_flux = self._tfa_model.solver.configuration.tolerances.feasibility

        sigma = sink.flux
        lump_dict = dict()

        if (sigma > epsilon_flux ):
            lump_dict[sink] = sigma
            for rxn in self._not_blocked_rxns:
                if (self._activation_vars[rxn].variable.primal < epsilon_int):
                    lump_dict[rxn] = rxn.flux
        else:
            self._tfa_model.logger.info('Metabolite {} is not able to produce'
                                        'by environment or external medium'.format(met_BBB.id))
        return lump_dict


    def _build_lump_sum_flux(self, met_BBB, sink):
        """
        This function uses the current solution of self._tfa_model

        :param met_BBB:
        :param sink:
        :return:
        """

        # import pdb; pdb.set_trace()
        epsilon_flux = self._tfa_model.solver.configuration.tolerances.feasibility
        # rxns=self._not_blocked_rxns+[sink]
        # fluxes=np.array([r.flux for r in rxns])
        # non_zero_fluxes=fluxes[abs(fluxes)>epsilon_flux]
        # active_rxns=rxns[abs(fluxes)>epsilon_flux]
        #
        # try:
        #     lump_dict = dict(zip(active_rxns, non_zero_fluxes))
        # except:
        #     import pdb; pdb.set_trace()

        sigma = sink.flux
        lump_dict = dict()
        if (sigma > epsilon_flux ):
            lump_dict[sink] = sigma
            for rxn in self._not_blocked_rxns:
                if (abs(rxn.flux)>epsilon_flux):
                    lump_dict[rxn] = rxn.flux
        else:
            self._tfa_model.logger.info('Metabolite {} is not able to produce'
                                        'by environment or external medium'.format(met_BBB.id))
        return lump_dict


    def _lump_min_plus_p(self, met_BBB, sink, p, force_solve):
        """

        :param met_BBB:
        :param sink:
        :param force_solve:
        :return:
        """

        epsilon = self._tfa_model.solver.configuration.tolerances.integrality

        try:
            max_lumps =self._param_dict['max_lumps_per_BBB']
        except KeyError:
            # TODO: Put a warning
            max_lumps=10

        lumps = list()

        with self._tfa_model as model:
            activation_vars = model.get_variables_of_type(FluxKO)

            # Solve a first time, obtain minimal subnet
            # import pdb; pdb.set_trace()
            model.slim_optimize()
            max_deactivated_rxns = model.objective.value

            # Add constraint forbidding subnets bigger than p
            expr = symbol_sum(activation_vars)

            # The lower bound is the max number of deactivated, minus p
            # Which allows activating the minimal number of reactions, plus p
            # lb = max_deactivated_rxns - p
            model.add_constraint(kind=ForbiddenProfile,
                                 hook = model,
                                 id_ = 'MAX_DEACT_{}'.format(met_BBB.id),
                                 expr = expr,
                                 lb=0,
                                 ub=len(model.reactions)
                                 # lb = lb,
                                 # ub = max_deactivated_rxns,
                                 )

            n_deactivated_reactions = max_deactivated_rxns

            # While loop, break on infeasibility
            while len(lumps)<max_lumps:

                try:
                    this_lump = self._lump_one_per_bbb(met_BBB, sink, force_solve)
                except (InfeasibleExcept, TimeoutExcept) as e:
                    if force_solve:
                        pass
                    elif len(lumps) == 0:
                        # No solution AND no lump found
                        raise e
                if model.solver.status != OPTIMAL:
                    break
                elif  isinstance(this_lump, list) and this_lump[1]=='sum_flux':
                    this_lump=this_lump[0]
                    lumps.append(this_lump)
                    return lumps
                elif this_lump is None:
                    # Since the solver is optimal, and we caught optim errors before,
                    # Then the BBB is simply produced in enough quantity by the core
                    break

                lumps.append(this_lump)

                # Add constraint forbidding the previous solution

                is_inactivated = [x for x in activation_vars
                               if abs(x.variable.primal-1) < 2*epsilon]

                expr = symbol_sum(is_inactivated)
                model.add_constraint(kind=ForbiddenProfile,
                                     hook = model,
                                     id_ = '{}_{}_{}'.format(met_BBB.id,
                                                             n_deactivated_reactions,
                                                             len(lumps)),
                                     expr = expr,
                                     lb = max_deactivated_rxns-p-self._param_dict['diverge_min'],
                                     ub = n_deactivated_reactions-self._param_dict['diverge_min'],
                                     )

        # TODO: Update of dynamic properties not handled yet
        # upon exiting context manager
        model.repair()
        return lumps

    def sum_reactions(rxn_dict, id_ = 'summed_reaction', epsilon = 1e-9):
        """
        Keys are reactions
        Values are their multiplicative coefficient
        """
        stoich = defaultdict(int)

        for rxn,flux in rxn_dict.items():
            for x, coeff in rxn.metabolites.items():
                stoich[x.id] += coeff * flux

        gpr = ') and ('.join(x.gene_reaction_rule for x in rxn_dict if x.gene_reaction_rule)

        gpr = ('(' + gpr + ')') if gpr else ''

        stoich = trim_epsilon_mets(stoich, epsilon=epsilon)

        new = Lump(id_ = id_,
                   metabolites = stoich,
                   subnetwork = {x.id:v for x,v in rxn_dict.items()},
                   gene_reaction_rule=gpr)

        return new



    def compute_mins_network(self, force_solve=False, method='OnePerBBB'):
        """
        For each BBB (reactant of the biomass reaction), add the corresponding sink to the model, then optimize and
        lump the result into one lumped reaction
        :param force_solve: Indicates whether the computations must continue when one lumping yields a status "infeasible"
        :return: The dict {BBB: lump} containing every lumped reactions, associated to their BBBs
        """

        # # Must be called before optimization
        # self._tfa_model.convert()
        # self._tfa_model.objective_direction = 'min'

        epsilon = self._tfa_model.solver.configuration.tolerances.feasibility

        # self._tfa_model_save=self._tfa_model
        the_method = disambiguate(method)
        print('Min network method detected: {}'.format(the_method))

        # dict: {metabolite: lumped_reaction}
        lumps = {}

        self._tfa_model.objective_direction = 'max'

        sink_iter = tqdm(self._sinks.items(), desc = 'met')

        for met_BBB, (sink_id, stoech_coeff) in sink_iter:

            # Cute stuff

            sink_iter.set_description('met={}'.format(met_BBB.id[:10]))
            sink_iter.refresh()

            sink = self._tfa_model.reactions.get_by_id(sink_id)
            # Activate reaction by setting its lower bound
            prev_lb = sink.lower_bound
            min_prod = self._growth_rate * stoech_coeff
            sink.lower_bound = min_prod - epsilon

            ## apply sum of flux
            print('sink_lower_bound',sink.lower_bound)
            self._apply_sum_of_flux()
            ## variability_analysis
            try:
                fva=pickle.load(open(os.getcwd()+'/'+sink_id+'fva.pickle','rb'))
            except:
                fva = variability_analysis(self._tfa_model, kind='reactions')
                pickle.dump(fva,open(os.getcwd()+'/'+sink_id+'fva.pickle','wb'))


            epsilon = self._tfa_model.solver.configuration.tolerances.feasibility
            blocked_df = fva[ (fva.max(axis=1).abs()<1*epsilon)
                    & (fva.min(axis=1).abs()<1*epsilon)]


            ## non_blocked_reactions
            not_blocked_rxns=[]
            for rxn in self._tfa_model.reactions:
                if rxn.id not in blocked_df.index.to_list():
                    not_blocked_rxns.append(rxn)
            self._not_blocked_rxns=not_blocked_rxns
            rids=set([ r.id for r in not_blocked_rxns])
            comm=list(set(self.orig_rid)&rids)
            add_var_rxns=[]
            for id in comm:
                add_var_rxns.append(self._tfa_model.reactions.get_by_id(id))

            # tmp_model=self._prune_FBA_model(self._tfa_model,blocked_df.index)

             # lumpgem binary variables to deactivate non-core reactions. The reaction is deactivated when the value of
            # the variable is 1
            self._activation_vars = {rxn: self._tfa_model.add_variable(kind=FluxKO,
                                                                       hook=rxn,
                                                                       lb=0,
                                                                       ub=1,
                                                                       queue=False)
                                     for rxn in add_var_rxns}

            print('Number of reaction:',len(add_var_rxns))
            self._generate_usage_constraints()
            self._generate_objective()
            # self._sinks=self._prepare_Mts()
            ### remove sum_of_flux bound
            self._FluxsumVar.lb=0
            self._FluxsumVar.ub=0

            if the_method == 'oneperbbb':
                this_lump = self._lump_one_per_bbb(met_BBB, sink, force_solve)
                lumped_reactions = [this_lump] if this_lump is not None else list()
            elif the_method.startswith('min+'):
                try:
                    p = int(the_method.replace('min+',''))
                except ValueError:
                    raise ValueError('Min+p method must have p as an integer')
                lumped_reactions = self._lump_min_plus_p(met_BBB, sink, p, force_solve)
            elif the_method.startswith('min'):
                lumped_reactions = self._lump_min_plus_p(met_BBB, sink, 0, force_solve)
            else:
                raise ValueError('Lumping method not recognized: {}. '
                                 'Valid methods are '
                                 'OnePerBBB, Min, Min+p, p natural integer'
                                 .format(the_method))


            if not lumped_reactions:
                continue

            lumps[met_BBB] = lumped_reactions
            # Deactivating reaction by setting both bounds to 0
            sink.lower_bound = prev_lb
            # sink.knock_out()

        self.lumps = lumps
        return lumps
