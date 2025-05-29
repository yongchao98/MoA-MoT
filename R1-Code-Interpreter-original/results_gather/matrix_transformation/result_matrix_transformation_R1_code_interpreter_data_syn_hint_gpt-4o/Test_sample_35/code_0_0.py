def rotate_matrix_90_clockwise(matrix):
    # Number of rows and columns in the original matrix
    m = len(matrix)
    n = len(matrix[0])
    
    # Create a new matrix for the rotated result
    rotated_matrix = [[0] * m for _ in range(n)]
    
    # Fill the rotated matrix
    for i in range(m):
        for j in range(n):
            rotated_matrix[j][m - 1 - i] = matrix[i][j]
    
    return rotated_matrix

# Original matrix
matrix = [
    [4, 3, 3],
    [3, 5, 1]
]

# Rotate the matrix
rotated = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated)