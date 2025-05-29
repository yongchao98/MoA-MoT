import numpy as np

# Define the input matrix
matrix = np.array([
    [1, 3, 9, 8, 9, 8, 2, 4],
    [5, 4, 1, 7, 3, 7, 1, 1],
    [0, 6, 4, 7, 2, 0, 8, 8],
    [4, 9, 7, 5, 4, 9, 1, 9],
    [0, 5, 9, 3, 3, 6, 2, 5],
    [0, 5, 1, 2, 1, 9, 1, 3],
    [8, 8, 9, 2, 9, 8, 0, 7],
    [3, 2, 0, 6, 0, 2, 3, 2]
])

# Format each element to two decimal places
formatted_matrix = np.round(matrix, 2)

# Print the formatted matrix
print(formatted_matrix)