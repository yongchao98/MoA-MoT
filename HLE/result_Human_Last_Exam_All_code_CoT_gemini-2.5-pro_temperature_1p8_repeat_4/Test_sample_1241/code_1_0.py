import numpy as np

# Given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# We need to solve the system of linear equations A*P = b
# where P = [P0, P1, P2, P3]^T are the steady-state probabilities.
# The equations are:
# l01*P0 - l10*P1 = 0
# l12*P1 - (l21+l23)*P2 = 0
# l23*P2 - l31*P3 = 0
# P0 + P1 + P2 + P3 = 1

# Set up the matrix A
A = np.array([
    [lambda_01, -lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

# Set up the vector b
b = np.array([0, 0, 0, 1])

# Solve the system A*P = b for P
try:
    P = np.linalg.solve(A, b)
    P0, P1, P2, P3 = P[0], P[1], P[2], P[3]

    # Calculate the required sum
    sum_P0_P1 = P0 + P1

    # Print the results
    print("The steady-state probabilities are:")
    print(f"P0(inf) = {P0}")
    print(f"P1(inf) = {P1}")
    print(f"P2(inf) = {P2}")
    print(f"P3(inf) = {P3}")
    print("\nThe sum of all probabilities (should be 1):")
    print(f"{P0} + {P1} + {P2} + {P3} = {P0 + P1 + P2 + P3}")

    print("\nThe required sum is P0(inf) + P1(inf):")
    print(f"{P0} + {P1} = {sum_P0_P1}")

except np.linalg.LinAlgError:
    print("The system of equations does not have a unique solution.")
