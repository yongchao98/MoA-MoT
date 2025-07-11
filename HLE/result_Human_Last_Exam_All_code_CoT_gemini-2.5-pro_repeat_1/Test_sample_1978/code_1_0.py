import numpy as np

# The index of the problem is the number of linearly independent boundary conditions.
# We can group the conditions by the variable they affect: x_1, x_2, and x_2024.
# The total index is the sum of the number of independent conditions for each variable.

# --- Conditions for x_1 ---
# Equation 1: x_1(T) - x_1(0) = 0  => -1*x_1(0) + 1*x_1(T) = 0
# The variables are (x_1(0), x_1(T)). The coefficient matrix is:
M1 = np.array([
    [-1, 1]
])
rank1 = np.linalg.matrix_rank(M1)

# --- Conditions for x_2 ---
# Equation 2: x_2(T) - x_2(0) = 0        => -1*x_2(0) + 1*x_2(T) = 0
# Equation 3: 5*x_2(T) - 5*x_2(0) = 0      => -5*x_2(0) + 5*x_2(T) = 0
# Equation 4: 100*x_2(T) - 100*x_2(0) = 0  => -100*x_2(0) + 100*x_2(T) = 0
# Equation 5: 1000*x_2(0) - 1000*x_2(T) = 0 => 1000*x_2(0) - 1000*x_2(T) = 0
# The variables are (x_2(0), x_2(T)). The coefficient matrix is:
M2 = np.array([
    [-1, 1],
    [-5, 5],
    [-100, 100],
    [1000, -1000]
])
rank2 = np.linalg.matrix_rank(M2)

# --- Conditions for x_2024 ---
# Equation 6: 100*x_2024(T) - 100*x_2024(0) = 0 => -100*x_2024(0) + 100*x_2024(T) = 0
# The variables are (x_2024(0), x_2024(T)). The coefficient matrix is:
M3 = np.array([
    [-100, 100]
])
rank3 = np.linalg.matrix_rank(M3)

# The total index is the sum of the ranks from each independent group of variables.
total_index = rank1 + rank2 + rank3

# Print the final equation showing how the total index is calculated
# from the number of independent equations in each group.
print(f"{rank1} + {rank2} + {rank3} = {total_index}")
