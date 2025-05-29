import numpy as np

# Define the matrix
matrix = np.array([
    [7, 9, 9],
    [4, 5, 7],
    [9, 4, 4],
    [7, 4, 5],
    [1, 8, 2]
])

# Remove every 3rd column (1-indexed), which is index 2 in 0-indexed
final_matrix = np.delete(matrix, 2, axis=1)

# Print the final matrix
print(final_matrix)