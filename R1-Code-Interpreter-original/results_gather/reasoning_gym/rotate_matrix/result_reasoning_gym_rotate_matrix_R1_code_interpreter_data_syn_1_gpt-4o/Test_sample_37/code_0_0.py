def rotate_90_clockwise(matrix):
    # For a 2x2 matrix, the rotation can be done by rearranging the elements
    return [
        [matrix[1][0], matrix[0][0]],
        [matrix[1][1], matrix[0][1]]
    ]

# Original matrix
matrix = [
    [1, 2],
    [1, 9]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)