# -*- coding: utf-8 -*-
# @author: Covariance (Pavel Martynov), 2022

from gibbs_energies.predictor import Predictor

def test_Al2O3():
    predictor = Predictor(
        initial_formula   = 'Al2O3',
        H                 = -3.442, # eV/atom
        path_to_structure = 'data/POSCAR.mp-1143_Al2O3',
        path_to_masses    = 'data/masses.json',
        path_to_chempots  = 'data/gels.json',
    )
    
    predictions = [predictor.dG(T) for T in (300, 600, 900, 1200, 1500, 1800)]
    
    true = [
        -3.2573151324,
        -3.0740994827,
        -2.9094292843,
        -2.7293955330,
        -2.5410811875,
        -2.3453501092,
    ]
    
    for actual, expected in zip(predictions, true):
        assert abs(actual - expected) < 1e-10

    
def test_dMgAl2O4():
    predictor = Predictor(
        initial_formula   = 'MgAl2O4',
        H                 = -3.404,
        path_to_structure = None,
        path_to_masses    = 'data/masses.json',
        path_to_chempots  = 'data/gels.json',
    )

    # Assuming V tabulated somewhere else
    predictions = [predictor.dG(T, 9.7) for T in (300, 600, 900, 1200, 1500, 1800)]
    
    true = [
        -3.2298861014,
        -3.0568000641,
        -2.9016493462,
        -2.7304384906,
        -2.5323029584,
        -2.3040658496,
    ]
    
    for actual, expected in zip(predictions, true):
        assert abs(actual - expected) < 1e-10
