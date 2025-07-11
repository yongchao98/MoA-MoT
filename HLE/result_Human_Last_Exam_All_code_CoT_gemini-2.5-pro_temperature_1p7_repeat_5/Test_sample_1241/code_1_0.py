import numpy as np

# Given lambda values
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> +inf), the derivatives P_i'(t) are 0.
# Let p_i = P_i(+inf). This gives us a system of linear equations.
# We will use the following equations to build the system Ap = b:
# 1) From P_0'(t)=0: -lambda_01*p0 + lambda_10*p1 = 0
# 2) From P_2'(t)=0: lambda_12*p1 - (lambda_21 + lambda_23)*p2 = 0
# 3) From P_3'(t)=0: lambda_23*p2 - lambda_31*p3 = 0
# 4) Normalization condition: p0 + p1 + p2 + p3 = 1
# Note: The equation for P_1'(t) becomes redundant with the corrected system.

# Coefficient matrix A
A = np.array([
    [-lambda_01, lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

# Constant vector b
b = np.array([0, 0, 0, 1])

# Solve the system of linear equations for p = [p0, p1, p2, p3]
try:
    p = np.linalg.solve(A, b)
    p0 = p[0]
    p1 = p[1]
    result = p0 + p1

    # Print the values for the final equation as requested
    print(f"P_0(inf) = {p0}")
    print(f"P_1(inf) = {p1}")
    print(f"The final equation is: P_0(inf) + P_1(inf) = {p0} + {p1}")
    print(f"Result: {result}")

except np.linalg.LinAlgError:
    print("The system of equations has no unique solution.")
