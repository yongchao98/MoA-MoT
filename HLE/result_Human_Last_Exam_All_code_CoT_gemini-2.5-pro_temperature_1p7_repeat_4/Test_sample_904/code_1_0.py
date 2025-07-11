import sympy as sp

# Define symbols for the variables involved
# gamma represents the surface tension
# r represents the radial position
gamma, r = sp.symbols('gamma r')

# Based on the derivation from the linearized Young-Laplace equation in cylindrical coordinates,
# the governing equation is:
# gamma * d^2(xi)/dr^2 + (gamma / r) * d(xi)/dr - Delta_P_electrostatic(r) = 0
#
# Comparing this to the form A(r)*d^2(xi)/dr^2 + B(r)*d(xi)/dr + C(r, xi) = 0,
# we can identify the coefficients A(r) and B(r).

A_r = gamma
B_r = gamma / r

# Print the results
print("Based on the derivation, the coefficients of the governing linear equation are:")
print(f"A(r) = {A_r}")
print(f"B(r) = {B_r}")
