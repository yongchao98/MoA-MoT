import math

# 1. Define constants based on the problem and necessary approximations.
# B: Spectral Radiance (given)
# L: Wavelength (given as 500 nm)
# C: Speed of Light (standard value)
# K: Boltzmann Constant. Approximated to 1.375e-23 (11/8e-23) to satisfy
#    the `frac` data type constraints of the described Wuxing architecture.
B_val = 9.9e16
L_val = 500e-9
C_val = 3e8
K_val = 1.375e-23
TWO = 2
FOUR = 4

# 2. Calculate the terms of the equation.
numerator = B_val * (L_val**FOUR)
denominator = TWO * C_val * K_val
temperature = numerator / denominator

# 3. Calculate the final answer in thousands of Kelvin, rounded.
final_answer = round(temperature / 1000)

# 4. Print the final equation with each number, as requested.
# The format `value:.3e` displays numbers in scientific notation.
print("The final equation with the numbers used:")
print(f"Temperature = ({B_val:.1e} * ({L_val:.0e})^{FOUR}) / ({TWO} * {C_val:.0e} * {K_val:.3e})")
print(f"Temperature = {numerator:.4e} / {denominator:.3e}")
print(f"Temperature = {temperature:.0f} K")
print("")
# 5. Print the final answer.
print("The final answer in a thousand Kelvin (rounded) is:")
print(int(final_answer))