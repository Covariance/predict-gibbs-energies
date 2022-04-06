# predict-gibbs-energies

Repository associated with [Physical descriptor for the Gibbs energy of inorganic crystalline solids and temperature-dependent materials chemistry](https://www.nature.com/articles/s41467-018-06682-4) article used to predict the temperature-dependent Gibbs energies of inorganic crystalline solids.

This repository is a rewritten copy of original repo.

## Data files

### table-S1.csv

Table of 440 compounds used for training and testing the SISSO-learned descriptor.

### masses.json

Dictionary of atomic masses (amu) in format `{element : mass}`.

### gels.json
  
Dictionary of experimental Gibbs energies (chemical potentials) for the elements in format `{temperature (K) : {element : G (eV/atom)}}`.

### POSCAR.mp-1143_Al2O3
  
Structure file for **Al2O3** from [Materials Project](https://materialsproject.org/) to demonstrate use of descriptor.

## Implementing descriptor

### predictor.py

Contains class for implementing SISSO-learned descriptor; see comments within and arXiv link for detailsю

### formula.py

Contains chemical formula parser; see comments within for usage examples.

## Changelog

- Formula standardizer extracted into separate class, now it is much more time-efficient (from 9x to 28x on different benchmarks) and able to process more complicated formulas;

- Tests provided to cover formula standartization.

- Predictor refactored and type hints provided.
