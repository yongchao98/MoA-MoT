from scipy.optimize import linprog
import numpy as np

# Coefficients for the objective function (we don't need to optimize, so use zeros)
c = np.zeros(9)

# Coefficients for the inequality constraints
A_eq = [
    [0, 1, 0, 0, 0, 0, -2, 0, 0],  # B = 2G
    [0, 0, 0, 0, 0, -2.5, 0, 0, 1],  # I = 2.5F
    [0, 1, 0, 0, 0, 0, 0, 0, 1],  # B + I = 86
    [1, 0, 0, 0, 0, 0, 1, 0, 0],  # A + G = 12
    [0, 0, 0, 0, 0, 1, 0, 0, -1],  # I > F (I - F > 0)
    [0, 0, 0, 0, 0, 1, 0, 0, 1],  # F + I = 112
    [0, 0, 1, 0, 0, 0, 0, 1, 0],  # C + H = 140
    [1, 0, 0, 0, 0, 0, 0, -1, 0]  # A - H = -11
]

b_eq = [0, 0, 86, 12, 0, 112, 140, -11]

# Bounds for each variable (they must be one of the given values)
bounds = [(3, 3), (6, 6), (9, 9), (20, 20), (32, 32), (36, 36), (80, 80), (120, 120), (121, 121)]

# Solve the linear programming problem
result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')

# Extract the solution
solution = result.x

# Print the solution in the required format
print(f"<<<{list(map(int, solution))}>>>")