def rotate_matrix_90_clockwise(matrix):
    # Number of rows in the original matrix
    n = len(matrix)
    # Number of columns in the original matrix
    m = len(matrix[0])
    
    # Create a new matrix with dimensions m x n
    rotated_matrix = [['' for _ in range(n)] for _ in range(m)]
    
    # Fill the rotated matrix with the appropriate elements
    for i in range(n):
        for j in range(m):
            rotated_matrix[j][n-1-i] = matrix[i][j]
    
    return rotated_matrix

# Original matrix
matrix = [
    ['4', '1', '3'],
    ['4', '7', '1']
]

# Rotate the matrix
rotated = rotate_matrix_90_clockwise(matrix)

# Format the output as requested
output = ','.join([','.join(row) for row in rotated])
print(f"<<<{output}>>>")