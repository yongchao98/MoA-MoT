import numpy as np

# Define the given transition rates (lambda values)
l01 = 0.019
l10 = 0.65
l12 = 0.4
l21 = 0.392
l23 = 0.008
l31 = 0.008

# Set up the coefficient matrix A and the constant vector b for the system A*P = b.
# The system of linear equations for the steady-state probabilities (P') is:
# -l01*P0 + l10*P1 + 0*P2 + 0*P3 = 0
# l12*P1 - (l21 + l23)*P2 + 0*P3 = 0
# 0*P0 + 0*P1 + l23*P2 - l31*P3 = 0
# P0 + P1 + P2 + P3 = 1

A = np.array([
    [-l01, l10, 0, 0],
    [0, l12, -(l21 + l23), 0],
    [0, 0, l23, -l31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Solve the system of linear equations for P = [P0, P1, P2, P3]
try:
    P = np.linalg.solve(A, b)
    P0 = P[0]
    P1 = P[1]
    
    # Calculate the required sum
    result = P0 + P1

    # Print the values and the final equation
    print(f"The steady-state probabilities are:")
    print(f"P0 = {P0:.6f}")
    print(f"P1 = {P1:.6f}")
    print(f"P2 = {P[2]:.6f}")
    print(f"P3 = {P[3]:.6f}")
    print("\nThe final sum is P0(+\u221e) + P1(+\u221e):")
    print(f"{P0:.6f} + {P1:.6f} = {result:.6f}")

except np.linalg.LinAlgError:
    print("The system of equations does not have a unique solution.")
