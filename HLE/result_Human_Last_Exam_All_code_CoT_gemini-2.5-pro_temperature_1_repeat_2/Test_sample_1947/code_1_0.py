import numpy as np

# Matrix A represents the coefficients of c_i for each graph
# Rows correspond to K3, C4, P4, K4, Lollipop graph
# Columns correspond to parameters e, k, p, sum(deg C 2), sum(deg C 3)
A = np.array([
    [3, 1, 0, 3, 0],
    [4, 0, 4, 4, 0],
    [3, 0, 1, 2, 0],
    [6, 4, 12, 12, 4],
    [4, 1, 2, 5, 1]
])

# Vector b represents the number of walks for each graph
b = np.array([6, 0, 0, 72, 6])

# Solve the linear system A*c = b
try:
    c = np.linalg.solve(A, b)
    # Print the coefficients
    print(f"{c[0]},{c[1]},{c[2]},{c[3]},{c[4]}")
except np.linalg.LinAlgError:
    print("The system of equations is singular and cannot be solved.")
