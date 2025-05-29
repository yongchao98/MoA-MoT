def rotate_matrix_90_clockwise(matrix):
    # Number of rows
    n = len(matrix)
    # Number of columns
    m = len(matrix[0])
    
    # Create a new matrix for the rotated result
    rotated_matrix = [[0] * n for _ in range(m)]
    
    # Perform the rotation
    for i in range(n):
        for j in range(m):
            rotated_matrix[j][n - 1 - i] = matrix[i][j]
    
    return rotated_matrix

# Given matrix
matrix = [
    [4, 2, 7],
    [7, 8, 6]
]

# Rotate the matrix
rotated = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated)