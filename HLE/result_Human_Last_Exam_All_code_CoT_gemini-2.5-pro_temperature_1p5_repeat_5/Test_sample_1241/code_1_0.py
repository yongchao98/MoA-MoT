import numpy as np

# Given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# We set up the system of linear equations A * p = b for the steady-state probabilities p = [p0, p1, p2, p3].
# We use the first, third, and fourth steady-state equations, plus the normalization equation.
# Equation 1: -lambda_01*p0 + lambda_10*p1 = 0
# Equation 3: lambda_12*p1 - (lambda_21 + lambda_23)*p2 = 0
# Equation 4: lambda_23*p2 - lambda_31*p3 = 0
# Normalization: p0 + p1 + p2 + p3 = 1

A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Solve the system of linear equations
try:
    p = np.linalg.solve(A, b)
    p0 = p[0]
    p1 = p[1]
    
    # Calculate the required sum
    result = p0 + p1
    
    # Print the equation with the found values
    print(f"P0(+inf) + P1(+inf) = {p0} + {p1} = {result}")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
