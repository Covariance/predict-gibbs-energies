# -*- coding: utf-8 -*-
# @author: Covariance (Pavel Martynov), 2022

from typing import Dict

# No type hints in sortedcontainers!
from sortedcontainers import SortedDict # type: ignore

from gibbs_energies.formula import StandardFormula


def table_test_parse(cases: Dict[str, Dict[str, int]]):
    for input, output in cases.items():
        assert StandardFormula(input).set == SortedDict(output)


def table_test_to_string(cases: Dict[str, str]):
    for input, output in cases.items():
        assert StandardFormula(input).string == output


def test_parse_simple():
    table_test_parse({
        'H2O'   : {'H' : 2, 'O' : 1},
        'H2 O'  : {'H' : 2, 'O' : 1},
        'H O H' : {'H' : 2, 'O' : 1},
        'HOH'   : {'H' : 2, 'O' : 1},
    })


def test_parse_nested():
    table_test_parse({
        '(H)2O'                             : {'H' : 2, 'O' : 1},
        '(H2)O'                             : {'H' : 2, 'O' : 1},
        '((H)H)O'                           : {'H' : 2, 'O' : 1},
        '(H2O)'                             : {'H' : 2, 'O' : 1},
        '(((H2))) O'                        : {'H' : 2, 'O' : 1},
        '(' * 100 + 'H' + ')' * 100 + '2 O' : {'H' : 2, 'O' : 1},
        
        '((((H)2)2)2)2' : {'H' : 16},
        '(((H)2O)2O)2O' : {'H' : 8, 'O' : 7},
    })

    
def test_parse_original_examples():
    table_test_parse({
        'O Ti2'             : {'O' : 1, 'Ti' : 2},
        'CaTiO3'            : {'Ca' : 1, 'O' : 3, 'Ti' : 1},
        'Al(OH)3'           : {'Al' : 1, 'H' : 3, 'O' : 3},
        '(Al10S)(OH2)3NNe2' : {'Al' : 10, 'H' : 6, 'N' : 1, 'Ne' : 2, 'O' : 3, 'S' : 1},
        
        'OAl(OH)3' : {'Al' : 1, 'H' : 3, 'O' : 4} # We're also able to parse this
    })


def test_to_string_original_examples():
    table_test_to_string({
        'O Ti2'             : 'O1Ti2',
        'CaTiO3'            : 'Ca1O3Ti1',
        'Al(OH)3'           : 'Al1H3O3',
        '(Al10S)(OH2)3NNe2' : 'Al10H6N1Ne2O3S1',
        'OAl(OH)3)'         : 'Al1H3O4',
    })
