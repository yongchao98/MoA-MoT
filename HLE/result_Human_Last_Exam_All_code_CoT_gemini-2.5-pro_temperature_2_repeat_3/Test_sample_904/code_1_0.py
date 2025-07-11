import sympy as sp

# Define the symbols used in the problem
r, gamma = sp.symbols('r gamma', real=True, positive=True)
xi = sp.Function('xi')(r)
E = sp.Function('E')(r)
epsilon_a, epsilon_b = sp.symbols('epsilon_a epsilon_b')

# Define the coefficients A(r) and B(r) based on the derivation
A_r = r
B_r = 1

# Define the C(r) term for completeness of the equation
C_r_xi = -r * (epsilon_b - epsilon_a) * E**2 / (2 * gamma)

# Construct the full differential equation
governing_equation = sp.Eq(A_r * xi.diff(r, 2) + B_r * xi.diff(r, 1) + C_r_xi, 0)

# Print the results
print("The governing linear equation for the interfacial shape xi(r) is:")
print(f"A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0\n")
print("Based on the derivation, the coefficients are:")
print(f"A(r) = {A_r}")
print(f"B(r) = {B_r}")
print("\nThe full equation is:")
# Using sp.pretty_print for better formatting of the equation
sp.pretty_print(governing_equation)

# Final answer to be extracted
print("\n<<<")
print(f"A(r) = {A_r}")
print(f"B(r) = {B_r}")
print(">>>")