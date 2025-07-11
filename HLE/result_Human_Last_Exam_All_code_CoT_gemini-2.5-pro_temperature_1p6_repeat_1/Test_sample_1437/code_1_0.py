import sympy
from sympy import Symbol, limit, log, pprint

# Step 1: Explain the physical formula.
# The Lindhard polarization function Pi_0 in the static limit (omega=0) is:
# Pi_0(k, 0) = -D(epsilon_F) * L(k)
# where D(epsilon_F) is the density of states at the Fermi energy, and L(k) is a dimensionless function.
# The question asks for the numerical value associated with the function at k=0.
# This corresponds to evaluating the dimensionless part, L(k), in the limit k -> 0.

# Step 2: Define the dimensionless function L(x) where x = k / (2*k_F).
# L(x) = 0.5 + ((1 - x^2) / (4*x)) * ln(|(1 + x) / (1 - x)|)
# We need to find the limit of L(x) as x -> 0.

# Define the symbolic variable x
x = Symbol('x')

# Define the dimensionless Lindhard function L(x).
# For |x| < 1, ln(|(1+x)/(1-x)|) simplifies to log((1+x)/(1-x)).
L_x = 0.5 + (1 - x**2) / (4 * x) * log((1 + x) / (1 - x))

print("The Lindhard polarization function is expressed as Pi_0 = -D(eps_F) * L(x),")
print("where D(eps_F) is the density of states at the Fermi level and x = k / (2*k_F).")
print("\nThe dimensionless function L(x) is:")
pprint(L_x, use_unicode=True)

# Step 3: Compute the limit of L(x) as x approaches 0.
result = limit(L_x, x, 0)

print("\nTo evaluate this at k=0 (which means x=0), we compute the limit:")
print(f"lim(x->0) L(x) = {result}")

# Step 4: Present the final equation with the computed numerical value.
# The problem asks to output the numbers in the final equation.
# The equation is Pi_0(k->0, omega=0) = -D(epsilon_F) * [result]
print("\nTherefore, the Lindhard polarization function at zero frequency and zero momentum is:")
print(f"Pi_0(k=0, w=0) = -D(epsilon_F) * {result}")
print("\nThe requested numerical value, which is the value of the dimensionless part of the function, is 1.")
