import numpy as np

# Define the input matrix
matrix = np.array([
    [2, 7, 8, 7, 0, 0, 8, 3, 4],
    [8, 4, 0, 5, 0, 6, 2, 4, 7],
    [5, 0, 9, 9, 5, 9, 0, 4, 6],
    [4, 2, 6, 7, 0, 8, 2, 7, 8],
    [0, 6, 4, 1, 9, 4, 9, 4, 5],
    [1, 5, 3, 7, 9, 3, 9, 5, 2],
    [5, 0, 3, 7, 2, 9, 7, 5, 5],
    [4, 4, 6, 5, 9, 9, 2, 8, 3],
    [4, 7, 2, 4, 8, 1, 9, 7, 9]
])

# Format each element to two decimal places
formatted_matrix = np.array([[f"{element:.2f}" for element in row] for row in matrix])

# Print the formatted matrix
print(formatted_matrix)