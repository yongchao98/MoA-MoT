import math

# Pandora's observational data
# B is the spectral radiance in W/m^2/sr/m
B = 9.9e16
# lmbda is the wavelength in meters
lmbda = 500e-9

# Physical constants
# c is the speed of light in m/s
c = 2.998e8
# k is the Boltzmann constant in J/K
k = 1.38e-23

# The integer constants from the equation
const_2 = 2
const_1000 = 1000

# Calculate the temperature in Kelvin using the Rayleigh-Jeans approximation:
# T = (B * λ^4) / (2 * c * k)
temperature_in_kelvin = (B * (lmbda**4)) / (const_2 * c * k)

# Divide by 1000 to get the value in "thousand Kelvin"
temperature_in_thousand_kelvin = temperature_in_kelvin / const_1000

# Round the result to the nearest integer
final_answer = round(temperature_in_thousand_kelvin)

# As requested, output the full equation with each number, and the final rounded answer.
# The "^" symbol is used here to denote exponentiation for readability.
print(f"round( (B * λ^4) / (2 * c * k) / 1000 )")
print(f"round( ({B} * ({lmbda})^4) / ({const_2} * {c} * {k}) / {const_1000} ) = {int(final_answer)}")