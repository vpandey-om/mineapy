
""" Minimum network enrichment analysis
.. moduleauthor:: vikash pandey
"""



from setuptools import setup, find_packages
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


version_tag = '0.0.4'

setup(name='mineapy',
      version=version_tag,
      author='vikash pandey',
      author_email='vikash.pandey@umu.se',
      url='https://github.com/vpandey-om/mineapy/',
      download_url='https://github.com/vpandey-om/mineapy/archive/v'+version_tag+'-alpha.tar.gz',
      install_requires=['cobra==0.20.0',
                        'bokeh==2.2.3',
                        'networkx==2.5',
                        'optlang==1.4.4',
                        'pytest==4.6.11',
                        'scipy==1.6.3',
                        'tqdm==4.56.0',
                        'pytfa==0.9.3',
                        'statsmodels==0.12.1',
                        'sympy==1.6.1'],
      packages = find_packages(),
      python_requires='>=3.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='mineapy, Minimum network enrichment analysis in Python',
      long_description=long_description,
      long_description_content_type='text/markdown',
      keywords=['mineapy','minea','thermodynamics','enrichment analysis'],

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
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
      ],
     )




