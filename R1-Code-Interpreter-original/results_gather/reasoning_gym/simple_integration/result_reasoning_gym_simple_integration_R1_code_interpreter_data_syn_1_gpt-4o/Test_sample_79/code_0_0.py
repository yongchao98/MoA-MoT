# Python code to simplify the second term
from sympy import symbols, simplify

X = symbols('X')
term = -231*X**7/35
simplified_term = simplify(term)
print(simplified_term)