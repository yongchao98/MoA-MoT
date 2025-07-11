import math

# 1. Define constants from the problem and physics.
B = 9.9e16      # Spectral radiance in W/m^2*sr*m
lmbda = 500e-9  # Wavelength in meters (500 nm)
h = 6.626e-34   # Planck's constant in J*s
c = 2.998e8     # Speed of light in m/s
k = 1.381e-23   # Boltzmann constant in J/K

# 2. As per the plan, the temperature T can be calculated with the formula:
# T ≈ (B*λ^4 / (2*c*k)) / (1 - h*c^2 / (B*λ^5))
# We will print this equation with all the numbers substituted.

print(f"T = (({B} * ({lmbda})**4) / (2 * {c} * {k})) / (1 - (({h} * ({c})**2) / ({B} * ({lmbda})**5)))")

# 3. Calculate the two main parts of the formula: the Rayleigh-Jeans approximation and the correction term.
T_approx_numerator = B * (lmbda**4)
T_approx_denominator = 2 * c * k
T_approx = T_approx_numerator / T_approx_denominator

C1_numerator = h * (c**2)
C1_denominator = B * (lmbda**5)
C1 = C1_numerator / C1_denominator

# 4. Calculate the final temperature using the derived formula.
T_actual = T_approx / (1 - C1)

# 5. Round the result to the nearest thousand Kelvin.
# The method is to add 500 and then perform integer division by 1000.
T_kilo_kelvin_rounded = int((T_actual + 500) / 1000)

# 6. Print the final rounded answer.
print(f"\nFinal Answer: {T_kilo_kelvin_rounded}")