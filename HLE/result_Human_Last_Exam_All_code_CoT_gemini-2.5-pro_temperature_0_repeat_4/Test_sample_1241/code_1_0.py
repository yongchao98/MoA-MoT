import numpy as np

# Step 1: Define the given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# Step 2: Set up the system of linear equations for the steady state (p_i' = 0).
# We will use the following equations:
# 1. lambda_01*p0 - lambda_10*p1 = 0
# 2. lambda_12*p1 - (lambda_21 + lambda_23)*p2 = 0
# 3. lambda_23*p2 - lambda_31*p3 = 0
# 4. p0 + p1 + p2 + p3 = 1
# This can be written in the matrix form A*p = b.

# Step 3: Define the matrix A and vector b
A = np.array([
    [lambda_01, -lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Step 4: Solve the system of equations for p = [p0, p1, p2, p3]
try:
    p = np.linalg.solve(A, b)
    p0 = p[0]
    p1 = p[1]
    
    # Step 5: Calculate the required sum P0(+inf) + P1(+inf)
    result = p0 + p1
    
    # Print the values of each probability and the final sum as an equation
    print(f"The steady-state probability P0(+inf) is: {p0}")
    print(f"The steady-state probability P1(+inf) is: {p1}")
    print("\nThe final equation is:")
    print(f"{p0} + {p1} = {result}")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
