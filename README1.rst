MiNEApy
=====
|PyPI|  |license|

This software can be used to perform minimal network enrichment analysis.

Vikash Pandey and Vassily Hatzimanikatis. Investigating the deregulation of metabolic tasks via Minimum Network Enrichment Analysis (MiNEA) as applied to nonalcoholic fatty liver disease using mouse and human omics data.DOI:https://doi.org/10.1371/journal.pcbi.1006760

Requirements
------------

You will need to have `Git-LFS <https://git-lfs.github.com/>`_ in order to properly download some binary files:

.. code:: bash

    git clone https://github.com/vpandey-om/mineapy.git /path/to/mineapy
    cd /path/to/mineapy
    git lfs install
    git lfs pull

**This module was developed in Python 3.7, and it is recommended to run Python 3.7
to run commercial solvers such as Gurobi and CPLEX.**


This module requires
`COBRApy <https://github.com/opencobra/cobrapy/>`_, `pyTFA <https://github.com/EPFL-LCSB/pytfa/>`_, as well as
`optlang <https://github.com/biosustain/optlang>`_ to work
properly. The installer should take care of that for you. You might also
want to install a dedicated solver. GLPK, CPLEX and Gurobi are
supported.


Setup
=====


You will need to install `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_ and create an environment using python 3.7:

.. code:: bash

    conda create -n py37 python=3.7
    conda activate py37

Using pip

.. code:: bash

    pip install -r requirements.txt
   or
   pip install mineapy==0.0.3

or from source

.. code:: bash

    conda create -n py37 python=3.7
    conda activate py37

    # inside the py37 environment we can download from github and install.
    git clone https://github.com/vpandey-om/mineapy.git /path/to/mineapy
    cd /path/to/mineapy
    pip install -e .


Quick start
===========
Two tutorial files detail thoroughly normal usages of the mineapy
package. They can be found at:

    ::

        mineapy
        └── tutorials
            ├── e_coli_core.py
            └── e_coli_gem.py


Usage
=====

First, create your COBRApy model. Make sure to define the additional
values required by pyTFA, as said in the "Models" page of the
documentation.

If you already have a Matlab model with thermodynamic data, you might
want to use ``pytfa.io.import_matlab_model``. Otherwise, have a look at
the `COBRApy
documentation <https://cobrapy.readthedocs.io/en/latest/io.html#MATLAB>`__,
then add the required properties.

If you're using a specific solver, don't forget to tell COBRApy about it
by setting the ``solver`` property of your model to the name of your
solver. See the `COBRApy
documentation <https://cobrapy.readthedocs.io/en/latest/solvers.html>`__
for more information about this.

Example script
--------------

Here is an example script :

.. code:: python

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


.. |PyPI| image:: https://img.shields.io/pypi/v/mineapy.svg
   :target: https://pypi.org/project/mineapy/

.. |license| image:: http://img.shields.io/badge/license-APACHE2-blue.svg
   :target: https://github.com/vpandey-om/mineapy/blob/main/LICENSE

License
========

The software in this repository is put under an APACHE-2.0 licensing scheme - please see the `LICENSE <https://github.com/EPFL-LCSB/pytfa/blob/master/LICENSE.txt>`_ file for more details.
