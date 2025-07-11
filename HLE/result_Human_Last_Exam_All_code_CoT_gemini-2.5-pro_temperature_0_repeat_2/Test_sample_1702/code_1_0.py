import math

# Pandora star temperature calculation

# Plan:
# 1. Use the Rayleigh-Jeans approximation of Planck's Law, as the specified
#    environment lacks the 'exp' function required for the full law.
#    The formula is: T = (Bλ * λ^4) / (2 * c * k)
# 2. Use constants approximated to a precision compatible with the described
#    'frac' data type.
# 3. The final answer is rounded to the nearest thousand Kelvin.

# Constants from the problem and physics
# Bλ: Spectral radiance in W/m^2·sr·m
B_lambda = 9.9e16
# λ: Wavelength in meters
lmbda = 500e-9
# c: Speed of light in m/s (approximated)
c = 3e8
# k: Boltzmann constant in J/K (approximated)
k = 1.4e-23
# A constant for the number 2
two = 2

# Perform the calculation based on the Rayleigh-Jeans Law
# T = (Bλ * λ^4) / (2 * c * k)
numerator = B_lambda * (lmbda**4)
denominator = two * c * k
temperature = numerator / denominator

# Round the result to the nearest thousand Kelvin
final_answer = int(round(temperature / 1000))

# Output the equation with all the numbers used in the calculation
print(f"Calculating temperature T using the equation: T = (Bλ * λ^4) / (2 * c * k)")
print(f"T = ({B_lambda} * ({lmbda})^4) / ({two} * {c} * {k})")
print(f"T = {numerator} / {denominator}")
print(f"T ≈ {int(round(temperature))} K")
print(f"\nFinal Answer (in thousand Kelvin, rounded): {final_answer}")