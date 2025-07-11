import numpy as np

# Given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# For the steady-state solution, we set the derivatives to zero.
# This results in a system of linear equations. We can represent this
# in matrix form A*P = b, where P = [P0, P1, P2, P3]^T.
# We will use three of the balance equations and the normalization condition.
#
# 1. -lambda_01*P0 + lambda_10*P1 = 0
# 2. lambda_12*P1 - (lambda_21 + lambda_23)*P2 = 0
# 3. lambda_23*P2 - lambda_31*P3 = 0
# 4. P0 + P1 + P2 + P3 = 1
#
# The matrix A and vector b are:
A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Solve the system of linear equations A*P = b
try:
    P = np.linalg.solve(A, b)
    P0 = P[0]
    P1 = P[1]
    
    # Calculate the required sum
    result = P0 + P1
    
    # Print the result in the requested format
    print("The final equation is:")
    print(f"P0(inf) + P1(inf) = P0 + P1")
    print(f"                    = {P0} + {P1}")
    print(f"                    = {result}")

except np.linalg.LinAlgError:
    print("The system of equations could not be solved due to a singularity.")
