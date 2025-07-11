import sympy

# Define the symbols used in the problem
# r represents the radial position
# gamma represents the surface tension
r, gamma = sympy.symbols('r gamma')

# From the derivation based on the linearized Young-Laplace equation,
# we identify the coefficients A(r) and B(r) for the governing differential equation.

# A(r) is the coefficient of the second derivative term d²ξ/dr²
A_r = gamma

# B(r) is the coefficient of the first derivative term dξ/dr
B_r = gamma / r

# Print the expressions for the coefficients
print(f"A(r) = {A_r}")
print(f"B(r) = {B_r}")