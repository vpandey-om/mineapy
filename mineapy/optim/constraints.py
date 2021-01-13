# -*- coding: utf-8 -*-
"""
.. module:: MiNEApy
   :platform: Unix, Mac
   :synopsis:Minimum network enrichment analysis

.. moduleauthor:: Vikash Pandey

Constraints declarations

"""
from pytfa.optim import ReactionConstraint, GenericConstraint, \
    ModelConstraint, GeneConstraint

class CatalyticConstraint(ReactionConstraint):
    """
    Class to represent a enzymatic constraint
    """

    prefix = 'CC_'


class SumFlux(GenericConstraint):
    """
    Class to represent a forbidden net flux directionality profile
    Looks like:
    F_rxn_1 + B_rxn_2 + ... + F_rxn_n <= n-1
    """

    def __init__(self, model, expr, id_, **kwargs):

        GenericConstraint.__init__(self,
                                   id_=id_,
                                   expr=expr,
                                   model=model,
                                   **kwargs)

    prefix = 'SumFlux_'
