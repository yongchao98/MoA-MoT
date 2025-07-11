import sympy

# Define the symbols
L = sympy.Symbol('L')
n = sympy.Symbol('n')

# The terms for P(n)
p2_numerator = 3*L**2 - 2*L + 2
p2_denominator = 24*n**2
P2 = p2_numerator / p2_denominator

p3_numerator = L**3 - 2*L**2 + 2*L
p3_denominator = 48*n**3
P3 = p3_numerator / p3_denominator

# The formula for P(n)
P_n = P2 + P3

# Pretty print the formula
# The instruction "output each number in the final equation!" suggests printing the components.
# The most direct representation is the sum of the two rational expressions.

formula_str = f"P(n) = ({sympy.pretty(p2_numerator)}) / ({sympy.pretty(p2_denominator)}) + ({sympy.pretty(p3_numerator)}) / ({sympy.pretty(p3_denominator)})"
final_formula = f"P(n) = (3*L**2 - 2*L + 2)/(24*n**2) + (L**3 - 2*L**2 + 2*L)/(48*n**3)"

# Combining the fractions to a single fraction
# P_n_combined = sympy.simplify(P_n)
# P(n) = (L**3*n + 4*L**2*n - 2*L**2 + 4*L + 4*n)/(48*n**3) # this is less clear

# The clearest way to represent the answer is the sum of fractions as derived.
print("The formula for P(n) is:")
print(f"P(n) = (3*L**2 - 2*L + 2) / (24*n**2) + (L**3 - 2*L**2 + 2*L) / (48*n**3)")