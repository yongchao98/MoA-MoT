import sympy

# Define symbols
L, n = sympy.symbols('L n')

# The coefficients of the terms in P(n)
p2_numerator = 3*L**2 + 2*L - 2
p2_denominator = 24
p2_term = p2_numerator / (p2_denominator * n**2)

p3_numerator = L**3 + 2*L**2 - 2*L
p3_denominator = 48
p3_term = p3_numerator / (p3_denominator * n**3)

# The formula for P(n)
P_n = p2_term + p3_term

# Print the formula for P(n)
# The prompt requires printing each number in the final equation.
# We will format the output string to be clear.
formula_str = f"P(n) = ({p2_numerator})/({p2_denominator}*n**2) + ({p3_numerator})/({p3_denominator}*n**3)"

print("The formula for P(n) where L = ln(n) is:")
print(formula_str)

# For a cleaner mathematical representation:
# print(f"P(n) = {sympy.pretty(P_n)}")