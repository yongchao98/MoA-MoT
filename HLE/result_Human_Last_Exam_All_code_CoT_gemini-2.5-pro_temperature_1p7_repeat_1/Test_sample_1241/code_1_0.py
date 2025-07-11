import numpy as np

# Step 1: Define the given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# Step 2: Solve the system of linear equations for the steady-state probabilities (pi_0, pi_1, pi_2, pi_3)
# The equations are:
# -lambda_01*pi_0 + lambda_10*pi_1 = 0                    => pi_0 = (lambda_10 / lambda_01) * pi_1
# lambda_12*pi_1 - (lambda_21 + lambda_23)*pi_2 = 0      => pi_2 = (lambda_12 / (lambda_21 + lambda_23)) * pi_1
# lambda_23*pi_2 - lambda_31*pi_3 = 0                    => pi_3 = (lambda_23 / lambda_31) * pi_2
# pi_0 + pi_1 + pi_2 + pi_3 = 1

# Express pi_0, pi_2, pi_3 in terms of pi_1
# Let's calculate the coefficients that relate them to pi_1
# pi_0 = coeff_01 * pi_1
coeff_01 = lambda_10 / lambda_01
# pi_2 = coeff_21 * pi_1
coeff_21 = lambda_12 / (lambda_21 + lambda_23)
# pi_3 = (lambda_23 / lambda_31) * pi_2 = (lambda_23 / lambda_31) * coeff_21 * pi_1
coeff_31 = (lambda_23 / lambda_31) * coeff_21

# Step 3: Use the normalization condition to find pi_1
# pi_1 * (coeff_01 + 1 + coeff_21 + coeff_31) = 1
sum_of_coeffs = coeff_01 + 1 + coeff_21 + coeff_31
pi_1 = 1 / sum_of_coeffs

# Step 4: Calculate pi_0
pi_0 = coeff_01 * pi_1

# Step 5: Calculate the final requested sum pi_0 + pi_1
result = pi_0 + pi_1

# Output the final equation with the calculated numbers
print(f"The required sum is P0(+infinity) + P1(+infinity).")
print(f"The calculated steady-state probabilities are P0 = {pi_0} and P1 = {pi_1}.")
print(f"The final equation is:")
print(f"{pi_0} + {pi_1} = {result}")