def rotate_matrix_90_clockwise(matrix):
    # Number of rows and columns
    n = len(matrix)
    m = len(matrix[0])
    
    # Create a new matrix for the rotated result
    rotated_matrix = [['' for _ in range(n)] for _ in range(m)]
    
    # Perform the rotation
    for i in range(n):
        for j in range(m):
            rotated_matrix[j][n-1-i] = matrix[i][j]
    
    return rotated_matrix

# Original matrix
matrix = [
    ['7', '4', '6'],
    ['6', '9', '4']
]

# Rotate the matrix
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)