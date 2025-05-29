import numpy as np

# Given output matrix
output_matrix = np.array([
    [0.3369140625, 0.33251953125, 0.5283203125, 0.251708984375],
    [0.5146484375, 0.7861328125, 0.52001953125, 0.46630859375],
    [0.39404296875, 0.61865234375, 0.135009765625, 0.58544921875],
    [0.68701171875, 0.60400390625, 0.43505859375, 0.619140625]
])

# Calculate matrix A
matrix_A = 2 * output_matrix

# Convert to list of lists for JSON serializability
matrix_A_list = matrix_A.tolist()

# Assume matrix B is a zero matrix
matrix_B_list = np.zeros_like(output_matrix).tolist()

# Combine A and B into kc_matrices
kc_matrices = [matrix_A_list, matrix_B_list]

print(kc_matrices)