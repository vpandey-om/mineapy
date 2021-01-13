pyTFA
PyPI Documentation Status Build Status Codecov Codacy branch grade license

Thermodynamics-based Flux Analysis, in Python. Paper : Pierre Salvy, Georgios Fengos, Meric Ataman, Thomas Pathier, Keng C Soh, Vassily Hatzimanikatis. "pyTFA and matTFA: a Python package and a Matlab toolbox for Thermodynamics-based Flux Analysis" Bioinformatics (2018), bty499, DOI: https://doi.org/10.1093/bioinformatics/bty499

Implements: Christopher S. Henry, Linda J. Broadbelt, and Vassily Hatzimanikatis. "Thermodynamics-based metabolic flux analysis." Biophysical journal 92.5 (2007): 1792-1805. DOI: https://doi.org/10.1529/biophysj.106.093138

Requirements
You will need to have Git-LFS in order to properly download some binary files:

git clone https://github.com/EPFL-LCSB/pytfa.git /path/to/pytfa
cd /path/to/pytfa
git lfs install
git lfs pull
This module was developed in Python 3.5, and it is recommended to run Python 3.5 to run commercial solvers such as Gurobi and CPLEX. Other Python versions (2.7, 3.4) should also work (see the CI builds)

This module requires COBRApy, as well as optlang to work properly. The installer should take care of that for you. You might also want to install a dedicated solver. GLPK, CPLEX and Gurobi are supported.

Container-based install
You might want to use this program inside of a container. The docker/ subfolder has all the necessary information and source files to set it up.

Setup
This step is not required if you're using the container, which bundles all this.

You can install this module with pip:

For Python 3, you might have to use pip3 instead of pip

pip3 install pytfa
or from source

git clone https://github.com/EPFL-LCSB/pytfa.git /path/to/pytfa
pip3 install -e /path/to/pytfa
Quick start
Three tutorial files detail thoroughly normal usages of the pytfa package. They can be found at:

pytfa
└── tutorials
    ├── figure_paper.py
    ├── tutorial_basics.py
    └── tutorial_sampling.py
More information can be found here.

Documentation
Documentation is hosted at Read the Docs

Alternatively you can also generate the docs locally.

Make sure sphinx is installed, and install as well the theme (this is already bundled with the container):

pip install sphinx sphinx-rtd-theme
You can then generate the documentation with this command:

cd work/pytfa/doc && make html
The resulting HTML files will be located in work/pytfa/doc/_build.

Testing the code
We recommend using the Docker container for testing the code, as it comes with everything bundled.

Install pytest if you don't already have it (pip install pytest, already included in the container), then start the tests with the pytest command.

Usage
First, create your COBRApy model. Make sure to define the additional values required by pyTFA, as said in the "Models" page of the documentation.

If you already have a Matlab model with thermodynamic data, you might want to use pytfa.io.import_matlab_model. Otherwise, have a look at the COBRApy documentation, then add the required properties.

If you're using a specific solver, don't forget to tell COBRApy about it by setting the solver property of your model to the name of your solver. See the COBRApy documentation for more information about this.

Thermodynamic database
You also need a thermodynamic database. Use thermoDBconverter.py if you have a thermodynamic database from Matlab you wish to import to Python.

Thermodynamic databases are stored in .thermodb files and can be easily loaded with pytfa.io.load_thermoDB.

Example script
Here is an example script :

import pytfa
from pytfa.io import import_matlab_model, load_thermoDB


cobra_model = import_matlab_model('../models/small_yeast.mat')

thermo_data = load_thermoDB('../data/thermo_data.thermodb')

mytfa = pytfa.ThermoModel(thermo_data, cobra_model)
mytfa.solver = 'optlang-cplex'

## TFA conversion
mytfa.prepare()
mytfa.convert()

## Info on the model
mytfa.print_info()

## Optimality
tfa_solution = mytfa.optimize()
License
The software in this repository is put under an APACHE-2.0 licensing scheme - please see the LICENSE file for more details.
