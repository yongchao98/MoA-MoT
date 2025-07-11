import math

# The problem asks for the value of 1/p_n for n = 1000.
n = 1000

# Based on the derivation, the formula for 1/p_n is 4 * cos^2(pi / (n + 2)).
# We will use this formula to compute the value.

# These are the numerical components of the final equation for n=1000.
constant_multiplier = 4
numerator_of_angle = math.pi
denominator_of_angle = n + 2
power = 2

# Calculate the final value
angle = numerator_of_angle / denominator_of_angle
cos_of_angle = math.cos(angle)
result = constant_multiplier * (cos_of_angle ** power)

# Output the final equation and the computed value, as requested.
print(f"The final equation for 1/p_{n} is: {constant_multiplier} * (cos(pi / ({n} + 2)))^{power}")
print(f"For n={n}, the calculation is: {constant_multiplier} * (cos({numerator_of_angle} / {denominator_of_angle}))^{power}")
print(f"The final calculated value of 1/p_1000 is: {result}")