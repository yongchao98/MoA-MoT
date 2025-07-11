import numpy as np

# These are the transition rates given in the problem
lambda_01 = 0.019
lambda_10 = 0.65
lambda_12 = 0.4
lambda_21 = 0.392
lambda_23 = 0.008
lambda_31 = 0.008

# At steady state (t -> +inf), the derivatives P_i'(t) are zero.
# Let p_i = P_i(+inf). This gives a system of linear equations.
# We will use the following four independent equations to solve for p_i:
# 1.  lambda_01*p0 - lambda_10*p1 = 0
# 2.  lambda_12*p1 - (lambda_21 + lambda_23)*p2 = 0
# 3.  lambda_23*p2 - lambda_31*p3 = 0
# 4.  p0 + p1 + p2 + p3 = 1
# This can be written in the matrix form A * p = b, where p = [p0, p1, p2, p3]^T.

# We define the matrix A based on the coefficients of p_i in the equations.
A = np.array([
    [lambda_01, -lambda_10, 0, 0],
    [0, lambda_12, -(lambda_21 + lambda_23), 0],
    [0, 0, lambda_23, -lambda_31],
    [1, 1, 1, 1]
])

# We define the vector b based on the constant terms in the equations.
b = np.array([0, 0, 0, 1])

# We solve the system A*p = b for the vector of probabilities p.
try:
    p = np.linalg.solve(A, b)
    p0 = p[0]
    p1 = p[1]

    # The final result is the sum of P0(+inf) and P1(+inf).
    result = p0 + p1

    print("The final equation is:")
    print(f"{p0} + {p1} = {result}")

except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
