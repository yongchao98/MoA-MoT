import math

# This program calculates the temperature of the star Pandora based on its spectral radiance,
# adhering to the computational constraints of the Wuxing architecture.

# Plan:
# 1. Use the Rayleigh-Jeans Law approximation of Planck's Law, as the `exp` function is not available.
#    The formula to find temperature T is: T = (B * λ^4) / (2 * c * k)
# 2. Define all constants ensuring they are representable by the `frac` data type,
#    where the numerator `n` must be a signed char (-128 to 127).
# 3. Perform the calculation and round the final result to the nearest thousand Kelvin.

# Step 1: Define the constants based on the problem and `frac` constraints.

# B = Spectral Radiance = 9.9e16 W/m^2srm
# This can be represented as frac{n=99, d=10, e=16}, which is valid.
B = 9.9e16

# λ = Wavelength = 500 nm = 5e-7 m
# This can be represented as frac{n=5, d=1, e=-7}, which is valid.
lambda_val = 5e-7

# c = Speed of light ≈ 3e8 m/s
# This can be represented as frac{n=3, d=1, e=8}, which is valid.
c = 3e8

# k = Boltzmann constant ≈ 1.38e-23 J/K
# The representation 138/100e-23 is invalid as n=138 is too large for a signed char.
# We use the close approximation 1.375e-23, which is 11/8 * 10^-23.
# This can be represented as frac{n=11, d=8, e=-23}, which is valid.
k = 1.375e-23

# The constant 2
two = 2

# Step 2: Calculate the temperature T using the formula.
# Each number used in the final equation is explicitly defined above.
lambda_4 = lambda_val ** 4
numerator = B * lambda_4
denominator = two * c * k
temperature = numerator / denominator

# Step 3: Convert to thousand Kelvin and round to the nearest integer.
temp_in_thousand_K = temperature / 1000
rounded_answer = int(temp_in_thousand_K + 0.5)

# Step 4: Print the final answer.
print(rounded_answer)