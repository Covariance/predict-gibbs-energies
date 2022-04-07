#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author: Covariance (Pavel Martynov), 2022

from typing import List, Optional

import click

from gibbs_energies.predictor import Predictor

@click.command()
@click.option('--formula', required=True, help='Chemical formula to calculate dG for', type=str)
@click.option('--H', required=True, help='Formation enthalpy at 0 or 298 K [eV/atom]', type=float)
@click.option('--structure', help='Path to DFT-optimized geometry file', type=str)
@click.option('--vol', help='Volume per atom as float', type=float)
@click.option('--masses', default='data/masses.json', help='Path to .json with masses', type=str)
@click.option('--chempots', default='data/gels.json', help='Path to .json with chempots', type=str)
@click.option('--temp', multiple=True, default=(300, 600, 900, 1200, 1500, 1800), help='Temperatures to calculate dG at')
def main(
    formula  : str,
    h        : float,
    structure: Optional[str],
    vol      : Optional[float],
    masses   : str,
    chempots : str,
    temp     : List[str],
    ):
    
    if structure is None and vol is None:
        raise ValueError('One of --structure and --vol must be set.')
    
    predictor = Predictor(formula, h, structure, masses, chempots)
    
    print(f'Approximating dG for {formula}:')
    for T in temp:
        print(f'\t T = {str(T).ljust(5)}: dG = {predictor.dG(T, vol):.10f} eV/atom')


if __name__ == '__main__':
    main()
