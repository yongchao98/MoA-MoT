import numpy as np

# Let the coefficients be c = [c1, c2, c3, c4, c5]
# We form a system of linear equations Ac = b

# Matrix A is formed by the parameters of 5 simple graphs
# Rows correspond to graphs: P3, K3, P4, K1,3, C4
# Columns correspond to parameters: e, k, p, sum(deg C 2), sum(deg C 3)
A = np.array([
    [2, 0, 0, 1, 0],  # G1 = P3
    [3, 1, 0, 3, 0],  # G2 = K3
    [3, 0, 1, 2, 0],  # G3 = P4
    [3, 0, 0, 3, 1],  # G4 = K_1,3
    [4, 0, 4, 4, 0]   # G5 = C4
])

# Vector b is formed by the number of closed tree-like walks of length 6 (N6) for each graph
b = np.array([
    0,   # N6 for P3
    0,   # N6 for K3
    6,   # N6 for P4
    12,  # N6 for K_1,3
    24   # N6 for C4
])

# Solve the system of equations
try:
    c = np.linalg.solve(A, b)
    # Convert coefficients to integers for cleaner output
    c_int = [int(round(x)) for x in c]

    # Print the coefficients
    print(','.join(map(str, c_int)))

except np.linalg.LinAlgError:
    print("The system of equations could not be solved (singular matrix).")
