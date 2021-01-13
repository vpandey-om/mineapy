# -*- coding: utf-8 -*-
"""
.. module::
   :platform: Unix, Windows
   :synopsis: flux balance models accounting for expression, thermodynamics, and resource allocation constraints

.. moduleauthor::

Variables declarations

"""

from pytfa.optim.variables import GenericVariable, BinaryVariable, \
    ReactionVariable, ModelVariable, GeneVariable, get_binary_type, \
    ForwardUseVariable,BackwardUseVariable



class FluxSumVar(ModelVariable):
    """
    Class to represent a product A*B when performing linearization of the
    model
    """
    prefix = 'FSV_'
