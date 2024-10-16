# -*- coding: utf-8 -*-
"""
.. module:: pytfa
   :platform: Unix, Windows
   :synopsis: Thermodynamics-based Flux Analysis

.. moduleauthor:: pyTFA team

Thermodynamic cobra_model class and methods definition


"""

import re
from copy import deepcopy
from math import log

import pandas as pd
from cobra import Model

from pytfa.core.model import LCSBModel
from pytfa.thermo import std
# from .metabolite import MetaboliteThermo
# from .reaction import calcDGtpt_rhs, calcDGR_cues, get_debye_huckel_b
# from .utils import (
#     check_reaction_balance,
#     check_transport_reaction,
#     find_transported_mets,
#     get_reaction_compartment,
# )
from pytfa.optim.constraints import (
    SimultaneousUse,
    NegativeDeltaG,
    BackwardDeltaGCoupling,
    ForwardDeltaGCoupling,
    BackwardDirectionCoupling,
    ForwardDirectionCoupling,
    ReactionConstraint,
    MetaboliteConstraint,
    DisplacementCoupling,
)
from pytfa.optim.variables import (
    ThermoDisplacement,
    DeltaGstd,
    DeltaG,
    ForwardUseVariable,
    BackwardUseVariable,
    LogConcentration,
    ReactionVariable,
    MetaboliteVariable,
)
from pytfa.utils import numerics
from pytfa.utils.logger import get_bistream_logger
#
BIGM = numerics.BIGM
BIGM_THERMO = numerics.BIGM_THERMO
BIGM_DG = numerics.BIGM_DG
BIGM_P = numerics.BIGM_P
EPSILON = numerics.EPSILON
MAX_STOICH = 10


class ThermoModel_WithoutInfo(LCSBModel, Model):
    """
    A class representing a cobra_model with thermodynamics information

    """

    def __init__(self, model=Model(), name=None):

        """
        :param float temperature: the temperature (K) at which to perform the calculations
        :param dict thermo_data: The thermodynamic database
        :type temperature: float
        """

        LCSBModel.__init__(self, model, name)

        self.logger = get_bistream_logger('ME model' + str(self.name))
        #self._init_thermo()

        self.logger.info('# Model initialized ')

    def _convert_reaction(self, rxn, verbose):
        """

        :param rxn:
        :param verbose:
        :return:
        """
        FU_rxn = self.add_variable(ForwardUseVariable, rxn)
        BU_rxn = self.add_variable(BackwardUseVariable, rxn)

        # create the prevent simultaneous use constraints
        # SU_rxn: FU_rxn + BU_rxn <= 1
        CLHS = FU_rxn + BU_rxn
        self.add_constraint(SimultaneousUse, rxn, CLHS, ub=1)

        # create constraints that control fluxes with their use variables
        # UF_rxn: F_rxn - M FU_rxn < 0
        F_rxn = rxn.forward_variable
        CLHS = F_rxn - FU_rxn * BIGM
        self.add_constraint(ForwardDirectionCoupling, rxn, CLHS, ub=0)

        # UR_rxn: R_rxn - M RU_rxn < 0
        R_rxn = rxn.reverse_variable
        CLHS = R_rxn - BU_rxn * BIGM
        self.add_constraint(BackwardDirectionCoupling, rxn, CLHS, ub=0)

    def convert_withoutInfo(
        self, verbose=True
    ):
        """ Converts a cobra_model into a tFBA ready cobra_model by adding the
        thermodynamic constraints required

        .. warning::
            This function requires you to have already called
            :func:`~.pytfa.ThermoModel.prepare`, otherwise it will raise an Exception !

        """

        self.logger.info("# Model conversion starting...")

        ###########################################
        # CONSTANTS & PARAMETERS for tFBA problem #
        ###########################################

        # value for the bigM in big M constraints such as:
        # UF_rxn: F_rxn - M*FU_rxn < 0
        bigM = BIGM
        # Check each reactions' bounds
        for reaction in self.reactions:
            if (
                reaction.lower_bound < -bigM - EPSILON
                or reaction.upper_bound > bigM + EPSILON
            ):
                raise Exception("flux bounds too wide or big M not big enough")
            if reaction.lower_bound < -bigM:
                reaction.lower_bound = -bigM
            if reaction.upper_bound > bigM:
                reaction.upper_bound = bigM

        ###################
        # INPUTS & CHECKS #
        ###################

        ## For each reaction...
        for rxn in self.reactions:
            self._convert_reaction(
                rxn, verbose
            )

        # CONSISTENCY CHECKS

        # Creating the objective
        if len(self.objective.variables) == 0:
            self.logger.warning("Objective not found")

        self.logger.info("# Model conversion done.")
        self.logger.info("# Updating cobra_model variables...")
        self.repair()
        self.logger.info("# cobra_model variables are up-to-date")

    def print_info(self, specific=False):
        """
        Print information and counts for the cobra_model
        :return:
        """
        if not specific:
            LCSBModel.print_info(self)

        n_metabolites = len(self.metabolites)
        n_reactions = len(self.reactions)
        n_metabolites_thermo = len(
            # [
            #     x
            #     for x in self.metabolites
            #     if hasattr(x, "thermo") and x.thermo["id"]
            # ]
            self._var_kinds[LogConcentration.__name__]
        )
        n_reactions_thermo = len(
            # [
            #     x
            #     for x in self.reactions
            #     if x.id is not None
            #     and hasattr(x, "thermo")
            #     and x.thermo["computed"]
            # ]
            self._cons_kinds[ForwardDeltaGCoupling.__name__]
        )

        info = pd.DataFrame(columns=["value"])
        info.loc["num metabolites(thermo)"] = n_metabolites_thermo
        info.loc["num reactions(thermo)"] = n_reactions_thermo
        info.loc["pct metabolites(thermo)"] = (
            n_metabolites_thermo / n_metabolites * 100
        )
        info.loc["pct reactions(thermo)"] = (
            n_reactions_thermo / n_reactions * 100
        )
        info.index.name = "key"

        print(info)

    def __deepcopy__(self, memo):
        """

        :param memo:
        :return:
        """

        return self.copy()

    def copy(self):

        from ..io.dict import model_from_dict, model_to_dict
        from ..optim.utils import copy_solver_configuration

        dictmodel = model_to_dict(self)
        new = model_from_dict(dictmodel)

        copy_solver_configuration(self, new)

        return new
