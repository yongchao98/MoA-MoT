import numpy as np

# Define the given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> +inf), the derivatives P_i'(t) are 0.
# Let p_i = P_i(+inf). We get a system of linear equations.
# We solve it by expressing p0, p2, p3 in terms of p1.

# From 0 = -lambda_01*p0 + lambda_10*p1  => p0 = (lambda_10 / lambda_01) * p1
coeff_p0 = lambda_10 / lambda_01

# From 0 = lambda_12*p1 - (lambda_21 + lambda_23)*p2 => p2 = (lambda_12 / (lambda_21 + lambda_23)) * p1
coeff_p2 = lambda_12 / (lambda_21 + lambda_23)

# From 0 = lambda_23*p2 - lambda_31*p3 => p3 = (lambda_23 / lambda_31) * p2
# Substituting the expression for p2: p3 = (lambda_23 / lambda_31) * coeff_p2 * p1
coeff_p3 = (lambda_23 / lambda_31) * coeff_p2

# Now use the normalization equation: p0 + p1 + p2 + p3 = 1
# (coeff_p0 * p1) + p1 + (coeff_p2 * p1) + (coeff_p3 * p1) = 1
# p1 * (coeff_p0 + 1 + coeff_p2 + coeff_p3) = 1
sum_of_coeffs = coeff_p0 + 1 + coeff_p2 + coeff_p3

# Solve for p1
p1 = 1 / sum_of_coeffs

# Solve for p0
p0 = coeff_p0 * p1

# The required sum is p0 + p1
result = p0 + p1

# Print the final result in the requested format
print("The values of the stationary probabilities are:")
print(f"p0 = P0(+inf) = {p0}")
print(f"p1 = P1(+inf) = {p1}")
p2 = coeff_p2 * p1
p3 = coeff_p3 * p1
print(f"p2 = P2(+inf) = {p2}")
print(f"p3 = P3(+inf) = {p3}")
print("\nThe final equation is:")
print(f"{p0} + {p1} = {result}")