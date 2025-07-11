import numpy as np

# Given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> +inf), P_i'(t) = 0.
# We get a system of linear equations for the stationary probabilities p_i.
# We choose three independent equations from the system and the normalization equation.
# -l_01*p0 + l_10*p1 = 0
# l_12*p1 - (l_21+l_23)*p2 = 0
# l_23*p2 - l_31*p3 = 0
# p0 + p1 + p2 + p3 = 1

# This can be written in matrix form Ax = b
# where x = [p0, p1, p2, p3]^T
A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Solve the system of linear equations for p = [p0, p1, p2, p3]
try:
    p = np.linalg.solve(A, b)
    p0, p1, p2, p3 = p[0], p[1], p[2], p[3]

    # Calculate the required sum
    result = p0 + p1
    
    # Print the results
    print("The system of linear equations is:")
    print(f"{-lambda_01}*p0 + {lambda_10}*p1 = 0")
    print(f"{lambda_12}*p1 - {lambda_21 + lambda_23}*p2 = 0")
    print(f"{lambda_23}*p2 - {lambda_31}*p3 = 0")
    print("p0 + p1 + p2 + p3 = 1\n")
    
    print("Solving this system gives the stationary probabilities:")
    print(f"P0(inf) = {p0:.7f}")
    print(f"P1(inf) = {p1:.7f}")
    print(f"P2(inf) = {p2:.7f}")
    print(f"P3(inf) = {p3:.7f}\n")

    print("The required sum is:")
    # Final output of the equation with the numbers plugged in
    print(f"P0(inf) + P1(inf) = {p0:.7f} + {p1:.7f} = {result:.7f}")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
