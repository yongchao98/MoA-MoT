import sympy

# Define the symbols used in the expressions
# r: radial position
# gamma: surface tension
r, gamma = sympy.symbols('r gamma')

# From the derivation, the coefficients A(r) and B(r) of the linear governing equation are:
# A(r) is the coefficient of the second derivative term (d^2ξ/dr^2)
A_r = gamma

# B(r) is the coefficient of the first derivative term (dξ/dr)
B_r = gamma / r

# Print the results
print(f"The coefficient A(r) is the coefficient of the second derivative term.")
print(f"A(r) = {A_r}")
print("\n")
print(f"The coefficient B(r) is the coefficient of the first derivative term.")
print(f"B(r) = {B_r}")
