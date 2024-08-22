from setuptools import setup, find_packages
from setuptools.command.install import install
import sys

class CustomInstallCommand(install):
    """Custom handler for the 'install' command."""
    
    def run(self):
        # Print the message about R requirements
        print("\n**************************************************")
        print("IMPORTANT: This package requires installation of R and download of the BiasedUrn package in R")
        print("Tested with R 4.2.2 and BiasedUrn 2.0.11")
        print("**************************************************\n")
        
        # Continue with the standard installation
        install.run(self)
        
setup(
    name='step_enrichment',
    version='0.1.0',
    author='Nicholas Markarian',
    author_email='nmarkari@usc.edu',
    description='enrichment on steps in Gene Ontology Causal Activity Models',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/nmarkari/gocam_enrichment',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'step_enrichment': ['data/*.dat'],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pandas',
        'rpy2',
        'scipy',
        'tqdm',
        #This package also requires installation of R. Download the R package BiasedUrn package: https://cran.r-project.org/web/packages/BiasedUrn/BiasedUrn.pdf. 
        #This code was tested with R version 4.2.2 and BiasedUrn_2.0.11.
    ],
    cmdclass={
        'install': CustomInstallCommand,
    },
)
