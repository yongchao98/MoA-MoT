import numpy as np

# Define the given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> +inf), the derivatives P_i'(t) are zero.
# This gives a system of linear equations for the steady-state probabilities p_i.
# 1) 0 = -lambda_01*p0 + lambda_10*p1
# 2) 0 = lambda_01*p0 - (lambda_10 + lambda_12)*p1 + lambda_21*p2 + lambda_31*p3
# 3) 0 = lambda_12*p1 - (lambda_21 + lambda_23)*p2
# 4) 0 = lambda_23*p2 - lambda_31*p3
# 5) 1 = p0 + p1 + p2 + p3

# We can solve this system by expressing p0, p2, and p3 in terms of p1.

# From equation (1): p0 = (lambda_10 / lambda_01) * p1
coeff_p0_p1 = lambda_10 / lambda_01

# From equation (3): p2 = (lambda_12 / (lambda_21 + lambda_23)) * p1
coeff_p2_p1 = lambda_12 / (lambda_21 + lambda_23)

# From equation (4): p3 = (lambda_23 / lambda_31) * p2
# Substitute the expression for p2: p3 = (lambda_23 / lambda_31) * coeff_p2_p1 * p1
coeff_p3_p2 = lambda_23 / lambda_31
coeff_p3_p1 = coeff_p3_p2 * coeff_p2_p1

# Now use the normalization condition (5):
# p0 + p1 + p2 + p3 = 1
# (coeff_p0_p1 * p1) + p1 + (coeff_p2_p1 * p1) + (coeff_p3_p1 * p1) = 1
# p1 * (coeff_p0_p1 + 1 + coeff_p2_p1 + coeff_p3_p1) = 1

# Calculate the total coefficient for p1
total_coeff = coeff_p0_p1 + 1 + coeff_p2_p1 + coeff_p3_p1

# Solve for p1
p1 = 1 / total_coeff

# Calculate p0
p0 = coeff_p0_p1 * p1

# The quantity to find is p0 + p1
result = p0 + p1

# Print the final equation with the calculated values
print("The final equation is P_0(+inf) + P_1(+inf) = Result")
print(f"{p0:.5f} + {p1:.5f} = {result:.5f}")