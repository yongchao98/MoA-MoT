def rotate_360(matrix):
    # A 360-degree rotation results in the same matrix
    return matrix

# Define the matrix
matrix = [
    [3, 4, 7, 8],
    [8, 2, 9, 8],
    [5, 8, 7, 7],
    [4, 1, 4, 6]
]

# Rotate the matrix by 360 degrees
rotated_matrix = rotate_360(matrix)

# Print the result
print(rotated_matrix)