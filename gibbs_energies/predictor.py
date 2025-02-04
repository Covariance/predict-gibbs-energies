# -*- coding: utf-8 -*-
# @author: Chris, Tue May 15 14:55:39 2018
# @author: Covariance (Pavel Martynov), 2022

import math
from typing import Any, Iterable, Optional, Dict
import json
from itertools import combinations

from pymatgen.core.structure import Structure

from gibbs_energies.formula import StandardFormula

class Predictor(object):
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
        self.formula = StandardFormula(initial_formula)
        self.path_to_structure = path_to_structure
        self.path_to_masses = path_to_masses
        self.path_to_chempots = path_to_chempots


    @property
    def standardize_formula(self) -> str:
        """
        Returns:
            Nicely formatted and alphabetized formula.
        
        Read more about conversion in gibbs_energies.formula.StandardFormula.
        """
        return self.formula.string


    @property
    def atom_names(self) -> Iterable[str]:
        """
        Returns:
            List of elements in formula in alphabetized order.
        """
        return self.formula.set.keys()


    @property
    def atom_nums(self) -> Iterable[int]:
        """
        Returns:
            List of numbers of elements in formula in alphabetized order.
        """
        return self.formula.set.values()


    @property
    def atom_total(self) -> int:
        """
        Returns:
            Number of atoms in formula unit.
        """
        return sum(self.atom_nums)


    @property
    def mass_d(self) -> Dict[str, float]:
        """
        Returns:
            Parsed masses json in the following format:
            {el : atomic mass [amu]}
        """
        with open(self.path_to_masses) as f:
            return json.load(f)


    @property
    def Gi_d(self) -> Any: # Currently Any because of complex form of chempots
        """
        Returns:
            Parsed chempots json in the following format:
            {temperature (str) [K] : {el (str) : G_el(T) (float) [eV/atom]}}
        """        
        with open(self.path_to_chempots) as f:
            return json.load(f)


    @property
    def m(self) -> float:
        """
        Returns:
            Reduced mass.
        """
        if len(self.formula.set) == 1:
            raise ValueError(f'Descriptor should not be applied to unary compounds: {self.formula.string}')

        for name in self.atom_names:
            if name not in self.mass_d.keys():
                raise ValueError(f'Unable to find mass of {name} in provided table')
        
        masses = map(self.mass_d.__getitem__, self.atom_names)
        
        total = .0
        for (mass1, num1), (mass2, num2) in combinations(zip(masses, self.atom_nums), 2):            
            total += (num1 + num2) * mass1 * mass2 / (mass1 + mass2)

        return total / ((len(self.formula.set) - 1) * self.atom_total)


    def V(self, vol_per_atom: Optional[float] = None) -> float:
        """
        Arguments:
            vol_per_atom - calculated atomic volume [A^3/atom] or None if reading in structure file
        Returns:
            calculated atomic volume [A^3/atom]
        """
        if self.path_to_structure is not None:
            struct = Structure.from_file(self.path_to_structure)
            return struct.volume / len(struct) 
        else:
            if vol_per_atom is None:
                raise ValueError('vol_per_atom cannot be None of path to structure is not provided')

            return vol_per_atom


    def Gd_sisso(self, T: float, vol_per_atom: Optional[float] = None) -> float:
        """
        Arguments:
            T - temperature [K]
            vol_per_atom - calculated atomic volume [A^3/atom] or None if reading in structure file
        Returns:
            G^delta as predicted by SISSO-learned descriptor (float) [eV/atom]
        """
        if len(self.formula.set) == 1:
            return .0
        else:
            m = self.m
            V = self.V(vol_per_atom)
            return (-2.48e-4 * math.log(V) - 8.94e-5 * m / V) * T + 0.181 * math.log(T) - 0.882


    def summed_Gi(self, T: float) -> float:
        """
        Arguments:
            T - temperature [K]
        Returns:
            sum of the stoichiometrically weighted chemical potentials of the elements at T (float) [eV/atom]
        """
        total = .0
        for el, num in zip(self.atom_names, self.atom_nums):
            
            if str(T) not in self.Gi_d or el not in self.Gi_d[str(T)]:
                raise ValueError(f'No entry in chempots for {el} with T={T}')
            
            Gi = self.Gi_d[str(T)][el]
            total += num * Gi
        return total


    def G(self, T: float, vol_per_atom: Optional[float] = None) -> float:
        """
        Args:
            T - temperature [K]
            vol_per_atom - calculated atomic volume [A^3/atom] or None if reading in structure file
        Returns:
            Absolute Gibbs energy at T using SISSO-learned descriptor for G^delta (float) [eV/atom]
        """
        if len(self.formula.set) == 1:
            # For now the type hint is ignored, because sortedcontainers is not typed
            return self.Gi_d[str(T)][self.atom_names[0]] # type: ignore
        else:
            return self.H + self.Gd_sisso(T, vol_per_atom)


    def dG(self, T: float, vol_per_atom: Optional[float] = None) -> float:
        """
        Arguments:
            T - temperature [K]
            vol_per_atom - calculated atomic volume [A^3/atom] or None if reading in structure file
        Returns:
            Gibbs formation energy at T using SISSO-learned descriptor for G^delta (float) [eV/atom]
        """
        if len(self.formula.set) == 1:
            return .0
        else:
            return ((self.H + self.Gd_sisso(T, vol_per_atom)) * 96.485 * self.atom_total - 96.485 * self.summed_Gi(T)) / (self.atom_total * 96.485)
