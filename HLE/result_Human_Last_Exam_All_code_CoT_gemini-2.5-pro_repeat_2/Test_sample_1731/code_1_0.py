import math
# Using scipy.constants for higher precision of physical constants
import scipy.constants as const

# --- Explanation ---
# This script finds the equilibrium values for mean energy density (u) and entropy density (s)
# for a photon gas (the Bose case of light quanta).
# The foundation for this calculation is the principle of maximum entropy, which is justified
# by large deviation theorems (like Sanov's theorem) in statistical mechanics. The procedure is:
# 1. Maximize the entropy for a Bose gas with no particle number conservation (chemical potential = 0)
#    subject to a fixed total energy. This yields the Bose-Einstein distribution.
# 2. Integrate the energy over all states using the Bose-Einstein distribution and the density of states
#    for photons to find the total energy U. This results in the Stefan-Boltzmann law for energy density u = U/V.
# 3. Use the thermodynamic relation S = (4/3) * U / T for a photon gas to find the entropy S and
#    the entropy density s = S/V.

# --- Physical Constants (in SI units) ---
k_B = const.k  # Boltzmann constant (J/K)
h = const.h    # Planck constant (J*s)
c = const.c    # Speed of light (m/s)
pi = math.pi

# --- Calculation of Coefficients ---

# The mean energy density u(T) is given by the formula:
# u(T) = (8 * pi^5 * k_B^4 / (15 * h^3 * c^3)) * T^4
# We calculate the coefficient A = (8 * pi^5 * k_B^4) / (15 * h^3 * c^3)
coeff_A = (8 * pi**5 * k_B**4) / (15 * h**3 * c**3)

# The entropy density s(T) is given by the formula:
# s(T) = (32 * pi^5 * k_B^4 / (45 * h^3 * c^3)) * T^3
# This coefficient is simply (4/3) * coeff_A
coeff_B = (4/3) * coeff_A

# --- Output the Final Equations ---
# The prompt requires outputting each number in the final equation. We present the
# symbolic formulas with their integer and constant components, followed by the
# final calculated form.

print("Equilibrium Mean Energy Density (u) for a Photon Gas:")
# Printing the numbers in the equation symbolically
print("Formula: u(T) = (8 * \u03c0^5 * k_B^4) / (15 * h^3 * c^3) * T^4")
# Printing the final equation with the calculated coefficient
print(f"Result:  u(T) = {coeff_A:.4e} * T^4  (in J/m^3)")
print("-" * 60)
print("Equilibrium Entropy Density (s) for a Photon Gas:")
# Printing the numbers in the equation symbolically
print("Formula: s(T) = (32 * \u03c0^5 * k_B^4) / (45 * h^3 * c^3) * T^3")
# Printing the final equation with the calculated coefficient
print(f"Result:  s(T) = {coeff_B:.4e} * T^3  (in J/(m^3*K))")