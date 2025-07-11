import numpy as np

# Given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# We express pi_0, pi_2, and pi_3 in terms of pi_1.
# pi_0 = c0 * pi_1
# pi_2 = c2 * pi_1
# pi_3 = c3_2 * pi_2 = c3_2 * c2 * pi_1 = c3 * pi_1
c0 = lambda_10 / lambda_01
c2 = lambda_12 / (lambda_21 + lambda_23)
c3_2 = lambda_23 / lambda_31
c3 = c3_2 * c2

# From the normalization equation: pi_0 + pi_1 + pi_2 + pi_3 = 1
# (c0 * pi_1) + pi_1 + (c2 * pi_1) + (c3 * pi_1) = 1
# pi_1 * (c0 + 1 + c2 + c3) = 1
sum_of_coeffs = c0 + 1 + c2 + c3
pi_1 = 1 / sum_of_coeffs

# Now calculate pi_0
pi_0 = c0 * pi_1

# The required value is pi_0 + pi_1
result = pi_0 + pi_1

# Print the final equation with the calculated values
# Using f-strings to format the numbers to a reasonable precision
print(f"The steady-state probabilities P0(+inf) and P1(+inf) are:")
print(f"P0(+inf) = {pi_0}")
print(f"P1(+inf) = {pi_1}")
print(f"The sum is P0(+inf) + P1(+inf) = {pi_0} + {pi_1} = {result}")
