
# MiNEApy

![PyPI](https://img.shields.io/pypi/v/mineapy.svg) ![license](http://img.shields.io/badge/license-APACHE2-blue.svg)

MiNEApy is a Python package for performing Minimum Network Enrichment Analysis (MiNEA), with a focus on enrichment analysis in biological systems.

Vikash Pandey and Vassily Hatzimanikatis. Investigating the deregulation of metabolic tasks via Minimum Network Enrichment Analysis (MiNEA) as applied to nonalcoholic fatty liver disease using mouse and human omics data. [DOI: https://doi.org/10.1371/journal.pcbi.1006760](https://doi.org/10.1371/journal.pcbi.1006760)

## Requirements

You will need to have [Git-LFS](https://git-lfs.github.com/) installed to properly download some binary files:

```bash
git clone https://github.com/vpandey-om/mineapy.git /path/to/mineapy
cd /path/to/mineapy
git lfs install
git lfs pull
```

**This module was developed in Python 3.7, and it is recommended to run Python 3.7 to ensure compatibility with commercial solvers such as Gurobi and CPLEX.**

MiNEApy depends on several Python packages including:

- [COBRApy](https://github.com/opencobra/cobrapy/)
- [pyTFA](https://github.com/EPFL-LCSB/pytfa/)
- [optlang](https://github.com/biosustain/optlang/)

The installation process should handle these dependencies for you. However, if you are using a dedicated solver, you may need to install it separately. Supported solvers include Gurobi and CPLEX.

## Setup

You will need to install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) and create an environment using Python 3.7:

```bash
conda create -n py37 python=3.7
conda activate py37
```

### Using pip

To install MiNEApy with pip, you can use the `requirements.txt` file:

```bash
pip install -r requirements.txt
```

Alternatively, you can install the specific version directly from PyPI:

```bash
pip install mineapy==0.0.4
```

### Installing from Source

If you prefer to install from the source:

```bash
conda create -n py37 python=3.7
conda activate py37

# Inside the py37 environment, download from GitHub and install.
git clone https://github.com/vpandey-om/mineapy.git /path/to/mineapy
cd /path/to/mineapy
pip install -e .
```

## Quick Start

With your environment activated, install Jupyter Notebook:

```
conda install -c conda-forge notebook
```

Install the ipykernel package to add your environment as a kernel in Jupyter:

```
conda install ipykernel
```

Add the environment to Jupyter as a kernel:

```
python -m ipykernel install --user --name=py37 --display-name "Python 3.7 (py37)"

```
# MiNEA Package Tutorials

I have created three tutorials to demonstrate how to use the MiNEA package. These tutorials are available in both Jupyter notebook (`.ipynb`) format and Python script (`.py`) format for ease of use. The tutorials guide users step-by-step through running MiNEA for different case studies.

### Available Tutorials:

1. **_E. coli_ Core Metabolic Model:**
   - Jupyter Notebook: `tutorial_1_e_coli_core.ipynb`
   - Python Script: `tutorial_1_e_coli_core.py`

2. **_E. coli_ Genome-Scale Model (GEM):**
   - Jupyter Notebook: `tutorial_2_ecoli_gem.ipynb`
   - Python Script: `tutorial_2_ecoli_gem.py`

3. **_Plasmodium berghei_ Model:**
   - Jupyter Notebook: `turorial_3_ipbe_blood_model.ipynb`
   - Python Script: `turorial_3_ipbe_blood_model.py`


## Usage

First, create your COBRApy model. Make sure to define the additional values required by pyTFA, as mentioned in the "Models" page of the documentation.

If you already have a Matlab model with thermodynamic data, you might want to use `pytfa.io.import_matlab_model`. Otherwise, refer to the [COBRApy documentation](https://cobrapy.readthedocs.io/en/latest/io.html#MATLAB) to add the required properties.

If you're using a specific solver, don't forget to configure COBRApy by setting the `solver` property of your model. See the [COBRApy solvers documentation](https://cobrapy.readthedocs.io/en/latest/solvers.html) for more information.

## License

The software in this repository is licensed under the Apache-2.0 License. Please see the [LICENSE](https://github.com/EPFL-LCSB/pytfa/blob/master/LICENSE.txt) file for more details.
