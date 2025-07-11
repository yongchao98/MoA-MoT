import math

# --- Physical Constants (SI units) ---
# Boltzmann constant (J/K)
k_B = 1.380649e-23
# Planck constant (J*s)
h = 6.62607015e-34
# Reduced Planck constant (J*s)
hbar = h / (2 * math.pi)
# Speed of light (m/s)
c = 299792458

# --- Input Parameters ---
# Temperature in Kelvin (e.g., surface of the sun is ~5778 K)
T = 5778
# Volume in cubic meters (e.g., 1 cubic meter)
V = 1.0

# --- Calculations ---

# 1. Equilibrium Mean Energy (U)
# This is derived by integrating the Planck distribution over the density of states.
# The result is the Stefan-Boltzmann law for the energy of black-body radiation.

hbar_c = hbar * c
k_T = k_B * T

U_numerator = math.pi**2 * V * k_T**4
U_denominator = 15 * hbar_c**3
U = U_numerator / U_denominator

print("--- Equilibrium Mean Energy (U) for a Photon Gas ---")
print("Derived from maximizing Bose-Einstein entropy, this gives the total energy in a volume V.")
print("\nFormula:")
print("U = (pi^2 * V * (k_B * T)^4) / (15 * (hbar * c)^3)")

print("\nCalculation with values:")
print(f"U = (pi^2 * {V:.1f} m^3 * ({k_B:.4g} J/K * {T} K)^4) / (15 * ({hbar:.4g} J*s * {c:.4g} m/s)^3)")
print(f"U = (pi^2 * {V:.1f} * ({k_T:.4g})^4) / (15 * ({hbar_c:.4g})^3)")
print(f"U = {U_numerator:.4g} / {U_denominator:.4g}")
print(f"U = {U:.4g} Joules")


# 2. Equilibrium Entropy (S)
# The entropy can be derived from the energy using the thermodynamic relation S = 4U / 3T.

S = (4 * U) / (3 * T)

print("\n\n--- Equilibrium Entropy (S) for a Photon Gas ---")
print("Derived from the thermodynamic relation for a photon gas.")
print("\nFormula:")
print("S = (4 * U) / (3 * T)")
print("S = (4 * pi^2 * V * k_B^4 * T^3) / (45 * (hbar * c)^3)")


print("\nCalculation with values:")
print(f"S = (4 * {U:.4g} J) / (3 * {T} K)")
print(f"S = {S:.4g} J/K")
