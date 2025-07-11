import numpy as np

# Given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> +inf), the derivatives P_i'(t) are 0.
# This gives a system of linear equations.
# We express P0, P2, and P3 in terms of P1.

# From P_3'(t) = 0 => lambda_23*P2 = lambda_31*P3
# P3 = (lambda_23 / lambda_31) * P2
# Since lambda_23 = lambda_31 = 0.008, P3 = P2.

# From P_2'(t) = 0 => lambda_12*P1 = (lambda_21 + lambda_23)*P2
# P2 = (lambda_12 / (lambda_21 + lambda_23)) * P1
# lambda_21 + lambda_23 = 0.392 + 0.008 = 0.4
# Since lambda_12 = 0.4, P2 = (0.4 / 0.4) * P1 = P1.
# Therefore, P1 = P2 = P3.

# From P_0'(t) = 0 => lambda_01*P0 = lambda_10*P1
# P0 = (lambda_10 / lambda_01) * P1
coeff_p0_p1 = lambda_10 / lambda_01

# The normalization condition is P0 + P1 + P2 + P3 = 1.
# Substitute the expressions in terms of P1:
# (coeff_p0_p1 * P1) + P1 + P1 + P1 = 1
# (coeff_p0_p1 + 3) * P1 = 1
P1 = 1 / (coeff_p0_p1 + 3)

# Now calculate P0
P0 = coeff_p0_p1 * P1

# The required sum is P0 + P1
total_prob = P0 + P1

# Print the final equation with the calculated values
print(f"The steady-state probabilities are P0 = {P0} and P1 = {P1}.")
print(f"The required sum is P0 + P1.")
print(f"{P0} + {P1} = {total_prob}")