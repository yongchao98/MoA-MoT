def rotate_matrix_720_degrees(matrix):
    # Since 720 degrees is equivalent to two full rotations (360 * 2),
    # the matrix remains unchanged.
    return matrix

# Define the input matrix
input_matrix = [
    [1, 3, 7, 4],
    [5, 2, 8, 2],
    [2, 0, 4, 2],
    [1, 0, 7, 7]
]

# Rotate the matrix by 720 degrees
output_matrix = rotate_matrix_720_degrees(input_matrix)

# Print the output matrix
print(output_matrix)