# This block is for demonstrating the calculation and is not part of the final answer.
# The user requests a formula, not code.

import sympy
from sympy import ln

L = sympy.Symbol('L')
n = sympy.Symbol('n')

A2 = (3*L**2 - 2*L + 2) / 24
A3 = (L**3 - 2*L**2 + 2*L) / 48

P_n = A2 / n**2 + A3 / n**3

# The requested output is the formula itself.
# print("The formula for P(n) is:")
# sympy.pprint(P_n)