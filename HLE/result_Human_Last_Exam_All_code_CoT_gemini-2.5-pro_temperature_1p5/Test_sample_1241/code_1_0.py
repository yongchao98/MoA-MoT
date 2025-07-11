import numpy as np

# Step 1: Define the given transition rates
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# Step 2: Set up the system of linear equations for the steady state (P_i' = 0)
# From P0'(t) = 0:  -lambda_01*P0 + lambda_10*P1 = 0
# From P2'(t) = 0:  lambda_12*P1 - (lambda_21 + lambda_23)*P2 = 0
# From P3'(t) = 0:  lambda_23*P2 - lambda_31*P3 = 0
# Normalization:    P0 + P1 + P2 + P3 = 1
# We can express this system in the matrix form A*P = b

# Step 3: Construct the matrix A and vector b
# P = [P0, P1, P2, P3]^T
A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Solve the system for the probability vector P
try:
    P_inf = np.linalg.solve(A, b)
    P0_inf = P_inf[0]
    P1_inf = P_inf[1]

    # Step 4: Calculate the required sum
    result = P0_inf + P1_inf

    # Print the individual probabilities and the final equation
    print(f"Calculated steady-state probabilities:")
    print(f"P0(+inf) = {P0_inf}")
    print(f"P1(+inf) = {P1_inf}")
    print(f"P2(+inf) = {P_inf[2]}")
    print(f"P3(+inf) = {P_inf[3]}")
    print("-" * 30)
    print("The final equation is:")
    print(f"P0(+inf) + P1(+inf) = {P0_inf} + {P1_inf}")
    print(f"Result = {result}")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
