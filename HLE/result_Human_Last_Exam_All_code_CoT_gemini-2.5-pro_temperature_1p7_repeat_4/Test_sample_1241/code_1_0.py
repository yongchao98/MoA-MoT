import numpy as np

# Given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# In the steady state (t -> +inf), P_i'(t) = 0. Let P_i = P_i(+inf).
# The system of differential equations becomes a system of linear equations:
# 1) 0 = -lambda_01 * P0 + lambda_10 * P1
# 2) 0 = lambda_01 * P0 - (lambda_10 + lambda_12) * P1 + lambda_21 * P2 + lambda_31 * P3
# 3) 0 = lambda_12 * P1 - (lambda_21 + lambda_23) * P2
# 4) 0 = lambda_23 * P2 - lambda_31 * P3
# 5) P0 + P1 + P2 + P3 = 1 (Normalization condition)

# From equation (1), we can express P0 in terms of P1:
# P0 = (lambda_10 / lambda_01) * P1
c0 = lambda_10 / lambda_01

# From equation (3), we can express P2 in terms of P1:
# P2 = (lambda_12 / (lambda_21 + lambda_23)) * P1
c2 = lambda_12 / (lambda_21 + lambda_23)

# From equation (4), we can express P3 in terms of P2:
# P3 = (lambda_23 / lambda_31) * P2
# Substituting the expression for P2:
# P3 = (lambda_23 / lambda_31) * c2 * P1
c3 = (lambda_23 / lambda_31) * c2

# Now substitute these into the normalization equation (5):
# c0*P1 + P1 + c2*P1 + c3*P1 = 1
# P1 * (c0 + 1 + c2 + c3) = 1
P1 = 1 / (c0 + 1 + c2 + c3)

# Now calculate P0 using the derived relation
P0 = c0 * P1

# The problem asks for the sum P0 + P1
sum_P0_P1 = P0 + P1

# Print the final equation with the calculated values
print(f"The stationary probabilities are calculated as:")
print(f"P0(inf) = {P0}")
print(f"P1(inf) = {P1}")
print(f"\nThe final sum is P0(inf) + P1(inf):")
print(f"{P0} + {P1} = {sum_P0_P1}")
