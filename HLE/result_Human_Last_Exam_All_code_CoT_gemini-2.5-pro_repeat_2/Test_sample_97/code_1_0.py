import math

# --- Given Constants ---
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

# --- Solving for [OH-] ---
# The charge balance equation is: 3[Al^3+] + [H+] = [Al(OH)4^-] + [OH^-]
# Substituting the equilibrium expressions gives:
# 3*K_sp/[OH-]^3 + K_w/[OH-] = K_sp*K_f*[OH-] + [OH-]
# Multiplying by [OH-]^3 gives:
# 3*K_sp + K_w*[OH-]^2 = K_sp*K_f*[OH-]^4 + [OH-]^4
# Rearranging into a quadratic equation for z = [OH-]^2:
# (1 + K_sp*K_f)*z^2 - K_w*z - 3*K_sp = 0

# Coefficients for the quadratic equation az^2 + bz + c = 0
a = 1 + K_sp * K_f
b = -K_w
c = -3 * K_sp

# --- Solve the quadratic equation for z = [OH-]^2 ---
# Using the quadratic formula: z = [-b +/- sqrt(b^2 - 4ac)] / 2a
discriminant = math.sqrt(b**2 - 4 * a * c)
# We take the positive root since z must be positive
z = (-b + discriminant) / (2 * a)

# --- Calculate concentrations ---
# [OH-] is the square root of z
OH_conc = math.sqrt(z)

# [Al^3+] = K_sp / [OH-]^3
Al3_conc = K_sp / (OH_conc**3)

# [Al(OH)4^-] = K_sp * K_f * [OH-]
AlOH4_conc = K_sp * K_f * OH_conc

# --- Calculate total solubility ---
# S = [Al^3+] + [Al(OH)4^-]
solubility = Al3_conc + AlOH4_conc

print("Step 1: Calculate the hydroxide ion concentration [OH-].")
print(f"[OH-] = {OH_conc:.2e} mol L^-1\n")

print("Step 2: Calculate the concentrations of the dissolved aluminum species.")
print(f"[Al^3+] = {Al3_conc:.2e} mol L^-1")
print(f"[Al(OH)4^-] = {AlOH4_conc:.2e} mol L^-1\n")

print("Step 3: Calculate the total solubility (S).")
print(f"S = [Al^3+] + [Al(OH)4^-]")
print(f"S = {Al3_conc:.2e} mol L^-1 + {AlOH4_conc:.2e} mol L^-1")
print(f"Total Solubility = {solubility:.2e} mol L^-1")
