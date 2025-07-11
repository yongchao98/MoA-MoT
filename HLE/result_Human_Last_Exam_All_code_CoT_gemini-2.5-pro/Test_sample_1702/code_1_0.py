import math

# Define the constants from the problem and physics
# B_lambda: spectral radiance from the problem
B_lambda = 9.9e16

# lambda_val: wavelength from the problem (500 nm converted to meters)
lambda_val = 500e-9

# two: the constant '2' from the equation
two = 2

# c: speed of light
c = 3e8

# k: Boltzmann constant
k = 1.38e-23

# Calculate the temperature in Kelvin using the Rayleigh-Jeans Law
# T = (B_lambda * lambda^4) / (2 * c * k)
temperature_kelvin = (B_lambda * (lambda_val**4)) / (two * c * k)

# Round the result to the nearest thousand Kelvin as requested
# For example, 51723 K -> round(51.723) -> 52
final_answer = int(round(temperature_kelvin / 1000))

# Per the instructions, output each number that is part of the final equation,
# followed by the final answer.
# The equation is T_rounded = round( (B_lambda * lambda_val^4) / (2 * c * k) / 1000 )
print(f"{B_lambda}")
print(f"{lambda_val}")
print(f"{two}")
print(f"{c}")
print(f"{k}")
print(f"{final_answer}")