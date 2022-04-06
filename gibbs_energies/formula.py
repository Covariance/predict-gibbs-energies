# -*- coding: utf-8 -*-
# @author: Covariance (Pavel Martynov), 2022

from sortedcontainers import SortedDict


class StandardFormula:
    """
    Class dedicated to parsing and storing chemical formulas. 
    
    Input formulas must follow the given grammar:
    
    formula := token | token WS formula
    token   := '(' formula ')' | token CNT | WORD
    
    WS   := arbitrary number of whitespaces, maybe zero
    WORD := name of a chemical element, namely capital letter
            followed by some number of lowercase letters, maybe zero
    CNT  := positive integer
    """
    def __init__(self, formula: str):
        self.initial_formula = formula
        
        self.set = StandardFormula._FormulaParser(formula).parse_formula()


    class _FormulaParser:
        """
        A single-pass chemical formula parser.
        """
        def __init__(self, source: str):
            self.source = source
            self.position = 0


        def current(self) -> str:
            """
            Returns next character from source,
            or '\0' if there's no characters left.
            """
            if self.position >= len(self.source):
                return '\0'
            
            return self.source[self.position]


        def next(self) -> str:
            """
            Returns current character from source,
            or '\0' if there's no characters left and
            moves pointer to the next one.
            """
            if self.position >= len(self.source):
                return '\0'
            
            self.position += 1
            
            return self.source[self.position - 1]


        def skip_ws(self) -> int:
            """
            Skips all whitespaces from the beginning
            of the source and returns the number
            of characters skipped.
            """
            skipped = 0
            
            while self.position < len(self.source) and self.source[self.position].isspace():
                self.position += 1
                skipped += 1
                
            return skipped


        def expect(self, expected: str):
            """
            Takes expected characters from source
            and throws ValueError on mismatch.
            """
            for ch in expected:
                if ch != self.next():
                    raise ValueError(f'Formula parse error: {ch} expected on position {self.position - 1}')


        def parse_formula(self) -> SortedDict:
            tokens = self.parse_token()
            
            self.skip_ws()
            
            while self.current() != '\0' and self.current() != ')':
                token = self.parse_token()
                
                for el, cnt in token.items():
                    if el not in tokens:
                        tokens[el] = cnt
                    else:
                        tokens[el] += cnt

                self.skip_ws()
            
            return tokens


        def parse_token(self) -> SortedDict:
            """
            Parses the token.
            """
            if self.current() == '(': # '(' token ')'
                self.next()
                token = self.parse_formula()
                self.expect(')')
            else: # WORD
                if not self.current().isupper():
                    raise ValueError(f'Expected capital letter at position {self.position}, found \'{self.current()}\'')
                
                word = self.next()
                while self.current().islower():
                    word += self.next()
                    
                token = SortedDict({word: 1})
              
            if self.current().isdigit(): # CNT is present
                digits = self.next()
                while self.current().isdigit():
                    digits += self.next()
                    
                try:
                    cnt = int(digits)
                except ValueError as err:
                    raise ValueError(f'Expected correct integer, but got \'{digits}\'') from err
            else: # CNT is not present
                cnt = 1
                
            
            for key in token.keys():
                token[key] = token[key] * cnt
            
            return token
    
    @property
    def string(self) -> str:
        """
        Returns formula formatted in a nice way, i.e.
        with elements alphabetically sorted and concatenated.
        """
        return ''.join(map(lambda x : x[0] + str(x[1]), self.set.items()))
