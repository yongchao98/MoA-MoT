import numpy as np

# Given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# We can solve this system using two methods.
# Method 1: Algebraic substitution (as planned above)

# Express P0, P2, P3 in terms of P1
# P0 = c0 * P1
# P2 = c2 * P1
# P3 = c3 * P1
c0 = lambda_10 / lambda_01
c2 = lambda_12 / (lambda_21 + lambda_23)
c3 = (lambda_23 / lambda_31) * c2

# The normalization equation is P0 + P1 + P2 + P3 = 1
# (c0*P1) + P1 + (c2*P1) + (c3*P1) = 1
# P1 * (c0 + 1 + c2 + c3) = 1
total_coeff = c0 + 1 + c2 + c3
P1 = 1 / total_coeff

# Calculate P0
P0 = c0 * P1

# Calculate the final sum
result = P0 + P1

# Print the final equation with the calculated values
print(f"The steady-state probabilities are P0(+inf) = {P0:.6f} and P1(+inf) = {P1:.6f}")
print(f"The required sum is P0(+inf) + P1(+inf):")
print(f"{P0:.6f} + {P1:.6f} = {result:.6f}")

# Method 2: Using numpy's linear algebra solver for verification
# The system of equations is:
# -l01*P0 + l10*P1 + 0*P2 + 0*P3 = 0
# 0*P0 + l12*P1 - (l21+l23)*P2 + 0*P3 = 0
# 0*P0 + 0*P1 + l23*P2 - l31*P3 = 0
# P0 + P1 + P2 + P3 = 1
# This can be written as A*P = b

A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Solve for P = [P0, P1, P2, P3]
P = np.linalg.solve(A, b)
P0_np, P1_np = P[0], P[1]
result_np = P0_np + P1_np

# print("\nVerification using numpy:")
# print(f"P0 = {P0_np:.6f}, P1 = {P1_np:.6f}")
# print(f"P0 + P1 = {result_np:.6f}")