def rotate_matrix_90_clockwise(matrix):
    # Number of rows and columns in the original matrix
    m, n = len(matrix), len(matrix[0])
    
    # Create a new matrix with dimensions n x m
    rotated_matrix = [[0] * m for _ in range(n)]
    
    # Fill the rotated matrix with the appropriate elements
    for i in range(m):
        for j in range(n):
            rotated_matrix[j][m - 1 - i] = matrix[i][j]
    
    return rotated_matrix

# Original matrix
matrix = [
    [8, 9, 1],
    [3, 8, 3]
]

# Rotate the matrix
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)