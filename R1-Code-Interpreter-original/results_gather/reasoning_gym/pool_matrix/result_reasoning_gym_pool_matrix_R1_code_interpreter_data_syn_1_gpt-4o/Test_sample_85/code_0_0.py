import numpy as np

# Define the input matrix
input_matrix = np.array([
    [0, 0, 1, 7],
    [3, 1, 5, 1],
    [0, 2, 9, 2],
    [6, 2, 5, 1],
    [1, 6, 7, 6]
])

# Format each element to two decimal places
output_matrix = np.array([[f"{element:.2f}" for element in row] for row in input_matrix])

print(output_matrix)