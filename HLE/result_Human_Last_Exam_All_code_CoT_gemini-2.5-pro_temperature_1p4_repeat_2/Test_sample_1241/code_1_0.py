import numpy as np

# Define the given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# We need to solve a system of linear equations for the steady-state probabilities p_i.
# The equations are:
# 1) -lambda_01*p0 + lambda_10*p1 = 0
# 2) lambda_12*p1 - (lambda_21 + lambda_23)*p2 = 0
# 3) lambda_23*p2 - lambda_31*p3 = 0
# 4) p0 + p1 + p2 + p3 = 1
# This can be written in matrix form A*p = b

# Set up the coefficient matrix A
A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

# Set up the result vector b
b = np.array([0, 0, 0, 1])

# Solve the system of linear equations A*p = b for p
try:
    p = np.linalg.solve(A, b)
    p0, p1, p2, p3 = p[0], p[1], p[2], p[3]
    
    # Calculate the required sum
    sum_p0_p1 = p0 + p1

    # Print the results
    print("The steady-state probabilities are:")
    print(f"P0(inf) = {p0:.6f}")
    print(f"P1(inf) = {p1:.6f}")
    print(f"P2(inf) = {p2:.6f}")
    print(f"P3(inf) = {p3:.6f}")
    print("\nThe problem is to find the value of P0(inf) + P1(inf).")
    print("\nFinal Equation:")
    print(f"P0(inf) + P1(inf) = {p0:.6f} + {p1:.6f} = {sum_p0_p1:.6f}")

except np.linalg.LinAlgError:
    print("The system of equations has no unique solution.")
