def rotate_matrix_360(matrix):
    # A 360-degree rotation results in the same matrix
    return matrix

# Define the input matrix
input_matrix = [
    [9, 8, 0],
    [5, 7, 8],
    [1, 7, 9]
]

# Rotate the matrix by 360 degrees
output_matrix = rotate_matrix_360(input_matrix)

# Print the output matrix
print(output_matrix)