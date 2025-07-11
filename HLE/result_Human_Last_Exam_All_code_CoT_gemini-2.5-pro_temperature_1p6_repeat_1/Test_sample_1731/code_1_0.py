import math

# --- Explanation ---
print("This script calculates the equilibrium mean energy (E) and entropy (S) for a photon gas in a defined volume and at a specific temperature.")
print("The calculations are based on the Stefan-Boltzmann law and thermodynamic relations derived from the principles of statistical mechanics, which are justified by large deviation theory.\n")
print("System Parameters:")
print("Volume (V) = 1.0 m^3")
print("Temperature (T) = 300.0 K\n")


# --- Constants ---
k_B = 1.380649e-23  # Boltzmann constant (J/K)
h = 6.62607015e-34   # Planck constant (J*s)
c = 299792458      # Speed of light (m/s)

# --- System Parameters ---
V = 1.0  # Volume (m^3)
T = 300.0  # Temperature (K)

# --- Calculation of Mean Energy (E) ---
print("1. Calculating the Mean Energy (E)")

# E = (8 * pi^5 * V * (k_B * T)^4) / (15 * h^3 * c^3)
E_numerator = 8 * (math.pi**5) * V * ((k_B * T)**4)
E_denominator = 15 * (h**3) * (c**3)
E = E_numerator / E_denominator

print("   Formula: E = (8 * \u03c0^5 * V * (k_B * T)^4) / (15 * h^3 * c^3)")
# Output each number in the final equation as requested
print(f"   E = (8 * {math.pi**5:.4f} * {V:.1f} * ({k_B:.6e} * {T:.1f})^4) / (15 * ({h:.6e})^3 * ({c:.6e})^3)")
print(f"   E = {E_numerator:.6e} / {E_denominator:.6e}")
print(f"   Equilibrium Mean Energy (E) = {E:.6e} Joules\n")

# --- Calculation of Entropy (S) ---
print("2. Calculating the Entropy (S)")

# S = (4 * E) / (3 * T)
S = (4 * E) / (3 * T)

print("   Formula: S = (4 * E) / (3 * T)")
# Output each number in the final equation as requested
print(f"   S = (4 * {E:.6e}) / (3 * {T:.1f})")
print(f"   Equilibrium Entropy (S) = {S:.6e} J/K")