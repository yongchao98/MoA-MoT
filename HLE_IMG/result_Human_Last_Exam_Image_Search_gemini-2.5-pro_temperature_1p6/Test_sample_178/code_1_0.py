import numpy as np

# Define the coefficient matrix A
A = np.array([
    [ 4, -1, -1,  0],
    [-1,  4,  0, -1],
    [-1,  0,  4, -1],
    [ 0, -1, -1,  4]
])

# Define the constant vector b
b = np.array([0.48, 1.46, -0.01, 0.98])

# Solve the system of linear equations A*T = b
T = np.linalg.solve(A, b)

# The result T will be in the order [T1, T2, T3, T4]
# print(T)