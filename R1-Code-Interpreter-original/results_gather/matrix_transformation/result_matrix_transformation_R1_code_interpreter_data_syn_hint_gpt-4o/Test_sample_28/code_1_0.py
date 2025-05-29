def rotate_matrix_90_clockwise(matrix):
    # Number of rows
    n = len(matrix)
    # Number of columns
    m = len(matrix[0])
    
    # Create a new matrix for the rotated result
    rotated_matrix = [['' for _ in range(n)] for _ in range(m)]
    
    # Perform the rotation
    for i in range(n):
        for j in range(m):
            rotated_matrix[j][n-1-i] = matrix[i][j]
    
    return rotated_matrix

# Given matrix
matrix = [
    ['6', '2', '8'],
    ['2', '9', '7']
]

# Rotate the matrix
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Flatten the rotated matrix and format the output
flattened_output = [item for sublist in rotated_matrix for item in sublist]
formatted_output = f"<<<{','.join(flattened_output)}>>>"

# Print the formatted output
print(formatted_output)