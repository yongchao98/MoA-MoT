import math

# Constants
K_sp = 5.3e-27
K_f = 1.1e31
Kw = 1.0e-14

# We solve for [OH-] by setting up and solving the charge balance equation,
# which results in a quadratic equation for y = [OH-]^2 of the form ay^2 + by + c = 0.
# (1 + K_sp*K_f)*[OH-]^4 - Kw*[OH-]^2 - 3*K_sp = 0

a = 1 + K_sp * K_f
b = -Kw
c = -3 * K_sp

# Solve the quadratic equation for y = [OH-]^2 using the quadratic formula.
# We take the positive root since concentration squared must be positive.
discriminant = b**2 - 4*a*c
y = (-b + math.sqrt(discriminant)) / (2*a)

# Calculate the hydroxide ion concentration
oh_conc = math.sqrt(y)

# Calculate the concentration of the aluminum-containing species
# [Al^3+] = K_sp / [OH-]^3
al3_conc = K_sp / (oh_conc**3)

# [Al(OH)4-] = (K_sp * K_f) * [OH-]
aloh4_conc = (K_sp * K_f) * oh_conc

# Total solubility (S) is the sum of the concentrations of all dissolved aluminum species
solubility = al3_conc + aloh4_conc

# Print the results
print("To find the solubility of Al(OH)3, we calculate the equilibrium concentrations of the dissolved aluminum species and sum them up.")
print("\nSolubility (S) = [Al^3+] + [Al(OH)4-]")
print(f"\nCalculated [Al^3+] = {al3_conc:.3e} mol L^-1")
print(f"Calculated [Al(OH)4-] = {aloh4_conc:.3e} mol L^-1")
print(f"\nS = {al3_conc:.3e} mol L^-1 + {aloh4_conc:.3e} mol L^-1")
print(f"\nThe total solubility of Al(OH)3 in pure water is {solubility:.3e} mol L^-1.")

# Final answer in the required format
# Rounding to two significant figures, consistent with the given constants.
final_answer = f"{solubility:.1e}"
print(f"\n<<<S = {final_answer} mol L^-1>>>")