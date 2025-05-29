def rotate_90_clockwise(matrix):
    # Number of rows
    n = len(matrix)
    # Number of columns
    m = len(matrix[0])
    
    # Create a new matrix for the rotated version
    rotated_matrix = [['' for _ in range(n)] for _ in range(m)]
    
    # Fill the rotated matrix
    for i in range(n):
        for j in range(m):
            rotated_matrix[j][n-1-i] = matrix[i][j]
    
    return rotated_matrix

# Original matrix
matrix = [
    ['%', '%', '^'],
    ['&', '%', '@', '@'],
    ['^', '&', '!', '@'],
    ['@', '$', '!', '!']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)