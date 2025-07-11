import numpy as np

# Define the transition rate coefficients
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> +inf), let p_i = P_i(+inf). The derivatives are zero.
# This gives a system of linear equations. We use three independent rate
# equations plus the normalization equation.
# 1. 0 = -lambda_01*p0 + lambda_10*p1
# 2. 0 = lambda_12*p1 - (lambda_21 + lambda_23)*p2
# 3. 0 = lambda_23*p2 - lambda_31*p3
# 4. 1 = p0 + p1 + p2 + p3
#
# We can write this system in matrix form A*p = b, where p = [p0, p1, p2, p3]^T.

A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Solve the system of linear equations for p
try:
    p = np.linalg.solve(A, b)
    p0 = p[0]
    p1 = p[1]

    # The value to be found is P0(+inf) + P1(+inf) = p0 + p1
    result = p0 + p1

    # Print the final equation with the calculated values
    print("The final equation is:")
    # The f-string formatting will display the full precision of the float variables.
    print(f"P0(inf) + P1(inf) = {p0} + {p1} = {result}")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
