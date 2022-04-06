# -*- coding: utf-8 -*-
"""
Created on Tue May 15 14:55:39 2018

@author: Chris
"""

from typing import Optional, Dict
import json
import re
from itertools import combinations
import math

from pymatgen.core.structure import Structure
import numpy as np

from .formula import StandardFormula

class PredictG(object):
    """
    Designed to take temperature-independent calculated or experimental data as input 
    and return the Gibbs formation energy at some temperature by applying
    a descriptor for the Gibbs energy of compounds
    """
    def __init__(self,
                 initial_formula  : str,
                 H                : float,
                 path_to_structure: Optional[str], 
                 path_to_masses   : str,
                 path_to_chempots : str
                 ):
        """
        Arguments:
        initial_formula - chemical formula (can be poorly formatted)
        H - formation enthalpy at 0 or 298 K [eV/atom]
        path_to_structure - path to DFT-optimized geometry file or None if providing volume per atom as float
        path_to_masses - path to .json with {el (str) : atomic mass (float) [amu]}
        path_to_chempots - path to .json with {temperature (str) [K] : {el (str) : G_el(T) (float) [eV/atom]}}
        """
        self.H = H
        self.initial_formula = initial_formula
        self.path_to_structure = path_to_structure
        self.path_to_masses = path_to_masses
        self.path_to_chempots = path_to_chempots
    
    @property
    def standardize_formula(self) -> str:
        """
        Returns nicely formatted and alphabetized formula.
        Read more about conversion in gibbs_energies.formula.StandardFormula.
        """
        return StandardFormula(self.initial_formula).string
        
    @property
    def atom_names(self):
        """
        Returns:
            list of alphabetized elements in formula (str)
        """
        formula = self.standardize_formula
        return re.findall('[A-Z][a-z]?', formula)
    
    @property
    def atom_nums(self):
        """
        Returns:
            list of alphabetized elements in formula (str)
        """    
        formula = self.standardize_formula
        return [int(num) for num in re.findall('\d+', formula)]
    
    @property
    def num_atoms(self):
        """
        Returns:
            number of atoms in formula unit (float)
        """
        return np.sum(self.atom_nums)

    @property
    def mass_d(self) -> Dict[str, float]:
        """
        Returns:
            {el (str) : atomic mass (float) [amu]} (dict)
        """
        with open(self.path_to_masses) as f:
            return json.load(f)
    
    @property
    def Gi_d(self):
        """
        Returns:
            {temperature (str) [K] : {el (str) : G_el(T) (float) [eV/atom]}} (dict)
        """        
        with open(self.path_to_chempots) as f:
            return json.load(f)
        
    @property
    def m(self):
        """
        Returns:
            reduced mass (float)
        """
        names = self.atom_names
        nums = self.atom_nums
        mass_d = self.mass_d
        num_els = len(names)
        num_atoms = np.sum(nums)        
        denom = (num_els - 1) * num_atoms
        if denom <= 0:
            print('descriptor should not be applied to unary compounds (elements)')
            return np.nan
        masses = [mass_d[el] for el in names]
        good_masses = [m for m in masses if not math.isnan(m)]
        if len(good_masses) != len(masses):
            for el in names:
                if math.isnan(mass_d[el]):
                    print('I dont have a mass for %s...' % el)
                    return np.nan
        else:
            pairs = list(combinations(names, 2))
            pair_red_lst = []
            for i in range(len(pairs)):
                first_elem = names.index(pairs[i][0])
                second_elem = names.index(pairs[i][1])
                pair_coeff = nums[first_elem] + nums[second_elem]
                pair_prod = masses[first_elem] * masses[second_elem]
                pair_sum = masses[first_elem] + masses[second_elem]
                pair_red = pair_coeff * pair_prod / pair_sum
                pair_red_lst.append(pair_red)
            return np.sum(pair_red_lst) / denom
            
    def V(self, vol_per_atom: Optional[float] = False):
        """
        Args:
            vol_per_atom (float or bool) - if reading in structure file, False; else the calculated atomic volume [A^3/atom]
        Returns:
            calculated atomic volume (float) [A^3/atom]
        """
        if self.path_to_structure != False:
            struct = Structure.from_file(self.path_to_structure )
            return struct.volume / len(struct) 
        else:
            return vol_per_atom
    
    def Gd_sisso(self, T, vol_per_atom=False):
        """
        Args:
            T (int) - temperature [K]
            vol_per_atom (float or bool) - if reading in structure file, False; else the calculated atomic volume [A^3/atom]        
        Returns:
            G^delta as predicted by SISSO-learned descriptor (float) [eV/atom]
        """
        if len(self.atom_names) == 1:
            return 0
        else:
            m = self.m
            V = self.V(vol_per_atom)
            return (-2.48e-4*np.log(V) - 8.94e-5*m/V)*T + 0.181*np.log(T) - 0.882
    
    def summed_Gi(self, T):
        """
        Args:
            T (int) - temperature [K]
        Returns:
            sum of the stoichiometrically weighted chemical potentials of the elements at T (float) [eV/atom]
        """
        names, nums = self.atom_names, self.atom_nums
        Gels = self.Gi_d
        els_sum = 0
        for i in range(len(names)):
            el = names[i]
            if el not in Gels[str(T)]:
                return np.nan
            num = nums[i]
            Gi = Gels[str(T)][el]
            els_sum += num*Gi
        return els_sum
    
    def G(self, T, vol_per_atom=False):
        """
        Args:
            T (int) - temperature [K]
            vol_per_atom (float or bool) - if reading in structure file, False; else the calculated atomic volume [A^3/atom]        
        Returns:
            Absolute Gibbs energy at T using SISSO-learned descriptor for G^delta (float) [eV/atom]
        """
        if len(self.atom_names) == 1:
            Gels = self.Gi_d
            el = self.atom_names[0]
            return Gels[str(T)][el]
        else:
            return self.H + self.Gd_sisso(T, vol_per_atom)
    
    def dG(self, T, vol_per_atom=False):
        """
        Args:
            T (int) - temperature [K]
            vol_per_atom (float or bool) - if reading in structure file, False; else the calculated atomic volume [A^3/atom]        
        Returns:
            Gibbs formation energy at T using SISSO-learned descriptor for G^delta (float) [eV/atom]
        """
        if len(self.atom_names) == 1:
            return 0.
        else:
            return ((self.H + self.Gd_sisso(T, vol_per_atom))*96.485*self.num_atoms - 96.485*self.summed_Gi(T)) / self.num_atoms / 96.485

def get_dGAl2O3_from_structure():
    """
    demonstration of how to get dG from optimized structure
    """
    print('------------------------------')    
    initial_formula = 'Al2O3'
    print('approximating dGf for %s...' % initial_formula)    
    path_to_structure = 'data/POSCAR.mp-1143_Al2O3'
    path_to_masses = 'data/masses.json'
    path_to_chempots = 'data/gels.json'
    H = -3.442 # eV/atom
    obj = PredictG(initial_formula,
                   H,
                   path_to_structure,
                   path_to_masses,
                   path_to_chempots)
    for T in [300, 600, 900, 1200, 1500, 1800]:
        print('T = %i K; dG = %.3f eV/atom' % (T, obj.dG(T=T, vol_per_atom=False)))
    print('------------------------------\n')
    return obj

def get_dMgAl2O4_without_structure():
    """
    demonstration of how to get dG from inputted volume per atom
    """
    print('------------------------------')
    initial_formula = 'MgAl2O4'
    print('approximating dGf for %s...' % initial_formula)
    path_to_structure = False
    V = 9.7 # A^3/atom (assuming tabulated somewhere, e.g. Materials Project)
    path_to_masses = 'data/masses.json'
    path_to_chempots = 'data/gels.json'
    H = -3.404 # eV/atom
    obj = PredictG(initial_formula,
                   H,
                   path_to_structure,
                   path_to_masses,
                   path_to_chempots)
    for T in [300, 600, 900, 1200, 1500, 1800]:
        print('T = %i K; dG = %.3f eV/atom' % (T, obj.dG(T=T, vol_per_atom=V)))
    print('------------------------------\n')
    return obj

def main():
    """
    run demonstrations
    Returns:
        PredictG objects
    """
    obj1 = get_dGAl2O3_from_structure()
    obj2 = get_dMgAl2O4_without_structure()
    return obj1, obj2

if __name__ == '__main__':
    obj1, obj2 = main()