def rotate_matrix_90_clockwise(matrix):
    # Transpose the matrix and then reverse each row
    return [list(reversed(col)) for col in zip(*matrix)]

# Original matrix
matrix = [
    ['4', '2', '9'],
    ['9', '4', '7']
]

# Rotate the matrix
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Flatten the rotated matrix to match the output format
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)