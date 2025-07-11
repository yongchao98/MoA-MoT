# The problem's constraints (no exp() function) require using the
# Rayleigh-Jeans approximation of Planck's Law: T ≈ (B * λ⁴) / (2 * c * k).
# This script calculates the temperature based on that formula.

# Given values
# Spectral Radiance B in W/m^2.sr.m
B = 9.9e16
# Wavelength lambda in meters (500 nm)
L = 500e-9

# Physical constants, simplified to values that would be representable
# in the described 'frac' type (e.g., with small integer numerators/denominators).
# Speed of light c in m/s
c = 3.0e8
# Boltzmann constant k in J/K (using 1.375e-23, which is 11/8 e-23)
k = 1.375e-23

# Calculate temperature T in Kelvin using the Rayleigh-Jeans formula.
# This is the "final equation" for the temperature calculation. The numbers
# used are B, L, c, and k, which are all defined as variables above.
temperature_K = (B * L**4) / (2 * c * k)

# The final step is to convert to thousand Kelvin and round.
# This equation uses the calculated temperature and the number 1000.
final_answer = round(temperature_K / 1000)

# Print the final result.
print(int(final_answer))