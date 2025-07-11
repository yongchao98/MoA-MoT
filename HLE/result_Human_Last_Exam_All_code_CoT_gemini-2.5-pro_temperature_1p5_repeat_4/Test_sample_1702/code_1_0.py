import math

# Define the constants based on the problem statement and known physical values.
# B_lambda: Spectral radiance in W/m^2/sr/m
# lmbda: Wavelength in meters (500 nm)
# c: Speed of light in m/s
# k: Boltzmann constant in J/K
B_lambda = 9.9e16
lmbda = 500e-9
c = 2.9979e8
k = 1.3806e-23

# The temperature T is calculated by rearranging the Rayleigh-Jeans Law:
# T = (B_lambda * lambda^4) / (2 * c * k)

# Calculate the numerator and denominator of the equation
numerator = B_lambda * (lmbda ** 4)
denominator = 2 * c * k

# Calculate the temperature in Kelvin
temperature_K = numerator / denominator

# The final answer must be in units of a thousand Kelvin, rounded to the nearest integer.
# For example, 51723 K -> round(51.723) -> 52.
final_answer = round(temperature_K / 1000)

# Print the final equation with its numerical values, as requested.
print("Using the Rayleigh-Jeans approximation: T = (Bλ * λ^4) / (2 * c * k)")
print(f"T = ({B_lambda:.1e} * ({lmbda:.1e})^4) / (2 * {c:.4e} * {k:.4e})")
print(f"Calculated Temperature: {int(round(temperature_K))} K")
print(f"Answer in thousand Kelvin (rounded): {final_answer}")