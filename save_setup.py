
""" Minimum network enrichment analysis
.. moduleauthor:: vikash pandey
"""

from setuptools import setup, find_packages
# import os
# from pip.req import parse_requirements
# from pip.download import PipSession

# __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
#
# def read_requirements():
#     '''parses requirements from requirements.txt'''
#     reqs_path = os.path.join(__location__, 'requirements.txt')
#     install_reqs = parse_requirements(reqs_path, session=PipSession())
#     reqs = [str(ir.req) for ir in install_reqs]
#     return reqs


version_tag = '1.1.0.0.0'


setup(
    name='mineapy',
    version=version_tag,
    author='vikash pandey',
    author_email='vikash.pandey@umu.se',
    url='https://github.com/vpandey-om/mineapy',
    packages=find_packages(include=['mineapy', 'mineapy.*']),
    install_requires=[
    'appdirs>=1.4.4',
    'attrs>=20.3.0',
    'bokeh>=2.2.3',
    'certifi>=2020.12.5',
    'chardet>=4.0.0',
    'cobra>=0.20.0',
    'colorama>=0.4.4',
    'commonmark>=0.9.1',
    'cplex>=12.10.0.0',
    'decorator>=4.4.2',
    'depinfo>=1.6.0',
    'diskcache>=5.1.0',
    'docloud>=1.0.375',
    'docplex>=2.19.202',
    'fastcache>=1.1.0',
    'future>=0.18.2',
    'gmpy2>=2.0.8',
     'h11>=0.11.0',
     'httpcore>=0.12.2',
    'httpx>=0.16.1',
    'idna>=2.10',
    'importlib-metadata>=3.3.0',
    'importlib-resources>=4.1.1',
    'iniconfig>=1.1.1',
    'Jinja2>=2.11.2'
    'MarkupSafe>=1.1.1',
    'mpmath>=1.1.0',
    'networkx>=2.5',
    'numpy>=1.19.4',
    'optlang>=1.4.4',
    'packaging>=20.8',
    'pandas>=1.2.0',
    'patsy>=0.5.1',
    'Pillow>=8.0.1',
    'pluggy>=0.13.1',
    'py>=1.10.0',
    'pydantic>=1.7.3',
    'Pygments>=2.7.3',
    'pyparsing>=2.4.7',
    'pytest>=6.2.1',
    'pytfa>=0.9.3',
    'python-dateutil>=2.8.1',
    'python-libsbml-experimental>=5.18.3',
    'pytz>=2020.5',
    'PyYAML>=5.3.1',
    'requests==2.25.1',
    #'rfc3986>=1.4.0'
    'rich>=6.2.0',
    'ruamel.yaml>=0.16.12',
    'ruamel.yaml.clib>=0.2.2',
    'scipy>=1.5.4',
    'six>=1.15.0',
    'sniffio>=1.2.0',
    'statsmodels>=0.12.1',
    #'swiglpk>=4.65.1',
    #'sympy @ file:///opt/concourse/worker/volumes/live/794b4585-6914-43d5-6dd6-9c84fb6d47ed/volume/sympy_1594236608820/work',
    'sympy==1.6.1',
    #'toml==0.10.2'
    #'tornado==6.1',
    'tqdm>=4.55.0',
    'typing-extensions>=3.7.4.3',
    'urllib3>=1.26.2',
    'zipp>=3.4.0'
    ],
    #python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
    python_requires='>=3.7',
    description='mineapy, Minimum network analysis in Python',

    license='Apache 2.0',

    # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
    'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Environment :: Console',

        # Pick your license as you wish (should match "license" above)
    'License :: OSI Approved :: Apache Software License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # 'Programming Language :: Python :: 3.5',
        # 'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    ]

)
