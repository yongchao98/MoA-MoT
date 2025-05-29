import numpy as np

# Define the input matrix
matrix = np.array([
    [6, 7, 5, 8, 8, 2, 1, 7, 2, 3],
    [7, 3, 9, 0, 1, 5, 9, 8, 4, 4],
    [8, 6, 3, 3, 7, 7, 1, 3, 5, 2],
    [2, 3, 7, 9, 3, 2, 2, 5, 0, 1],
    [7, 9, 3, 8, 5, 9, 3, 2, 6, 9],
    [9, 3, 5, 4, 5, 5, 8, 0, 1, 1],
    [4, 8, 9, 2, 5, 7, 8, 6, 4, 2],
    [9, 7, 3, 9, 0, 5, 5, 5, 0, 3],
    [2, 4, 9, 7, 6, 7, 7, 3, 4, 9]
])

# Format each element to two decimal places
formatted_matrix = np.array([[f"{element:.2f}" for element in row] for row in matrix])

# Print the formatted matrix
print(formatted_matrix)