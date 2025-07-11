import math

# This program calculates the temperature of the star Pandora based on its spectral radiance.
# Given the computational constraints of the Wuxing architecture (no exp or ln functions),
# the Rayleigh-Jeans approximation of Planck's Law is used.
# The formula is T = (B * λ^4) / (2 * c * k).

# --- Constants based on the problem statement ---

# Spectral radiance B = 9.9e16 W/m^2*sr*m
B = 9.9e16
B_str = "9.9e16"

# Wavelength λ = 500 nm = 500e-9 m
L = 500e-9
L_str = "500e-9"

# Speed of light c ≈ 3e8 m/s
c = 3.0e8
c_str = "3e8"

# Boltzmann constant k ≈ 1.4e-23 J/K.
# This value is used as it's a close approximation that fits the architecture's
# 'frac' type constraints (e.g., representable as 14/1e-24).
k = 1.4e-23
k_str = "1.4e-23"

# --- Calculation Steps ---

# Calculate T based on the rearranged Rayleigh-Jeans Law
L4 = L**4
numerator = B * L4
denominator = 2 * c * k
temperature = numerator / denominator

# The final answer must be in thousand Kelvin, rounded to the nearest integer.
final_answer = round(temperature / 1000)

# --- Output ---
# As requested, the following prints the components of the equation leading to the result.

print(f"Calculating Temperature (T) using Rayleigh-Jeans Law: T = (B * λ^4) / (2 * c * k)")
print(f"Using the following values:")
print(f"B (Spectral Radiance) = {B_str} W/m^2srm")
print(f"λ (Wavelength) = {L_str} m")
print(f"c (Speed of Light) = {c_str} m/s")
print(f"k (Boltzmann Constant) = {k_str} J/K")
print(f"")
print(f"Equation with numerical values:")
print(f"T = ({B_str} * ({L_str})^4) / (2 * {c_str} * {k_str})")
print(f"T = {numerator:.4e} / {denominator:.4e}")
print(f"T = {temperature:.2f} K")
print(f"")
print(f"Final Answer: Temperature in thousand Kelvin (rounded)")
print(f"Result = round({temperature:.2f} / 1000)")
print(f"Result = {final_answer}")
<<<737>>>