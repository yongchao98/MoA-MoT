import numpy as np

# Physical constants in SI units
k_B = 1.380649e-23  # Boltzmann constant (J/K)
h = 6.62607015e-34   # Planck constant (J*s)
c = 2.99792458e8    # Speed of light (m/s)
PI = np.pi

# --- Theory and Derivation Summary ---
# The equilibrium state of a system corresponds to the most probable macrostate.
# In large deviation theory (via the Boltzmann-Sanov theorem), this is the state
# whose empirical distribution minimizes the rate function (relative entropy).
# In physics, this is achieved by maximizing the Boltzmann entropy S = k * ln(W),
# where W is the number of microstates, subject to physical constraints.
#
# For a gas of photons (bosons with zero chemical potential), maximizing the entropy
# subject to a fixed total energy leads to the Bose-Einstein distribution.
#
# By integrating this distribution over the continuous density of states for photons
# in a volume V, we can find the total equilibrium energy U and entropy S.
# These derivations result in the famous Stefan-Boltzmann law for energy and a
# related law for entropy.

# 1. Calculate the radiation constant 'a' for the energy density equation u = a * T^4
# The constant 'a' is also known as the Stefan-Boltzmann constant divided by c/4.
# Formula: a = (8 * pi^5 * k_B^4) / (15 * h^3 * c^3)
a_rad_const = (8 * PI**5 * k_B**4) / (15 * h**3 * c**3)

# 2. Define a temperature for a worked example.
# We use the temperature of the Cosmic Microwave Background (CMB), a perfect photon gas.
T = 2.725  # Temperature in Kelvin

# 3. Calculate equilibrium values for energy density (U/V) and entropy density (S/V).
# Mean Energy Density (Stefan-Boltzmann Law): u = U/V = a * T^4
mean_energy_density = a_rad_const * T**4

# Equilibrium Entropy Density: s = S/V = (4/3) * a * T^3
entropy_density_coeff = (4 / 3) * a_rad_const
entropy_density = entropy_density_coeff * T**3

# 4. Print the results clearly, showing the final equations with their numerical coefficients.
print("--- Equilibrium Mean Energy for Photon Gas ---\n")
print("The mean energy density (U/V) is given by the Stefan-Boltzmann Law for radiation.")
print(f"Equation: U/V = a * T^4")
print(f"where the radiation constant a = (8 * pi^5 * k^4) / (15 * h^3 * c^3)")
print(f"Calculated value for a: {a_rad_const:.6e} J/(m^3 * K^4)\n")
# The request was to "output each number in the final equation"
print("Final equation with numerical value:")
print(f"U/V = ({a_rad_const:.6e}) * T^4  [J/m^3]\n")
print(f"For a temperature T = {T} K (Cosmic Microwave Background):")
print(f"Mean Energy Density = {mean_energy_density:.6e} J/m^3")
print("-" * 50)

print("\n--- Equilibrium Entropy for Photon Gas ---\n")
print("The equilibrium entropy density (S/V) is derived from thermodynamics (S = (U-F)/T).")
print(f"Equation: S/V = (4/3) * a * T^3")
print(f"The coefficient (4/3)*a is: {entropy_density_coeff:.6e} J/(m^3 * K^4)\n")
# Output each number in the final equation
print("Final equation with numerical value:")
print(f"S/V = ({entropy_density_coeff:.6e}) * T^3  [J/(K*m^3)]\n")
print(f"For a temperature T = {T} K (Cosmic Microwave Background):")
print(f"Entropy Density = {entropy_density:.6e} J/(K*m^3)")