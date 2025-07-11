import numpy as np

# Step 1: Define the transition rates from the problem description.
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# Step 2: Set up the system of linear equations for the steady state.
# At steady state (t -> +inf), P'(t) = 0. Let P_i(+inf) = p_i.
# We select three independent equations from the rate equations and add the normalization condition.
# The system in matrix form A*p = b is:
#
#  -lambda_01*p0 + lambda_10*p1 + 0*p2 + 0*p3 = 0
#  0*p0 + lambda_12*p1 - (lambda_21 + lambda_23)*p2 + 0*p3 = 0
#  0*p0 + 0*p1 + lambda_23*p2 - lambda_31*p3 = 0
#  1*p0 + 1*p1 + 1*p2 + 1*p3 = 1

# Step 3: Construct the coefficient matrix A and the constant vector b.
A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

b = np.array([0, 0, 0, 1])

# Step 4: Solve the system of equations for p = [p0, p1, p2, p3].
try:
    p = np.linalg.solve(A, b)
    p0, p1, p2, p3 = p[0], p[1], p[2], p[3]

    # Step 5: Calculate the required sum P0(+inf) + P1(+inf).
    result = p0 + p1

    # Step 6: Print the results, including the final equation with its components.
    print("The steady-state probabilities are:")
    print(f"P0(+inf) = {p0:.6f}")
    print(f"P1(+inf) = {p1:.6f}")
    print(f"P2(+inf) = {p2:.6f}")
    print(f"P3(+inf) = {p3:.6f}")
    print(f"\nCheck sum of probabilities: {np.sum(p):.6f}")
    
    print("\nThe final calculation is:")
    print(f"P0(+inf) + P1(+inf) = {p0:.6f} + {p1:.6f} = {result:.6f}")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
