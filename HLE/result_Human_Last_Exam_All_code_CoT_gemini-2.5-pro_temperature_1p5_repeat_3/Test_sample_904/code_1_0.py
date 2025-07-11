import sympy

# Define the symbols used in the problem
# r represents the radial position
# gamma represents the surface tension
r, gamma = sympy.symbols('r gamma')

# Based on the derivation, the coefficients A(r) and B(r) of the
# differential equation are determined by the linearized Young-Laplace equation
# in cylindrical coordinates.

# A(r) is the coefficient of the second derivative term d^2(xi)/dr^2
A_r = gamma

# B(r) is the coefficient of the first derivative term d(xi)/dr
B_r = gamma / r

# Print the results
print("The governing linear equation for the interfacial shape xi(r) has the form:")
print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0")
print("\nFrom the physical derivation based on surface tension, the coefficients are:")
print(f"A(r) = {A_r}")
print(f"B(r) = {B_r}")
