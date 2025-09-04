import numpy as np

# Define physical constants
h = 6.62607015e-34  # Planck's constant in J·s
c = 2.99792458e8    # Speed of light in m/s
k = 1.380649e-23     # Boltzmann constant in J/K

# Define the given parameters from the problem
T_nospots = 6000.0  # Temperature without spots in Kelvin
T_spots = 5500.0    # Temperature with spots in Kelvin
wavelength = 1448e-10 # Wavelength in meters (1448 Å)

# Step 1: Calculate the energy difference (ΔE) between the levels
delta_E = (h * c) / wavelength

# Step 2 & 3: Calculate the factor F.
# The ratio of populations is R = (g2/g1) * exp(-ΔE / (k*T)).
# The factor F is the ratio of the population ratios: F = R_nospots / R_spots.
# F = exp(-ΔE / (k*T_nospots)) / exp(-ΔE / (k*T_spots))
# This simplifies to F = exp( (ΔE / k) * (1/T_spots - 1/T_nospots) )
factor = np.exp((delta_E / k) * (1.0/T_spots - 1.0/T_nospots))

print(f"The energy difference (ΔE) is: {delta_E:.4e} J")
print(f"The factor by which the population ratio changes is: {factor:.2f}")