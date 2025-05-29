def rotate_matrix_90_clockwise(matrix):
    # Number of rows and columns in the original matrix
    m, n = len(matrix), len(matrix[0])
    
    # Create a new matrix with dimensions n x m
    rotated_matrix = [[0] * m for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(m):
        for j in range(n):
            rotated_matrix[j][m - 1 - i] = matrix[i][j]
    
    return rotated_matrix

# Original matrix
matrix = [
    [2, 1, 3],
    [3, 7, 7]
]

# Rotate the matrix
rotated = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated)