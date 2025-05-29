def rotate_matrix_90_counterclockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse the order of the rows
    rotated_matrix = transposed_matrix[::-1]
    return [list(row) for row in rotated_matrix]

# Original matrix
matrix = [
    [7, 8, 4, 7],
    [7, 6, 3, 7],
    [5, 9, 6, 3],
    [3, 0, 7, 8]
]

# Rotate the matrix 90 degrees counterclockwise
rotated_matrix = rotate_matrix_90_counterclockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)