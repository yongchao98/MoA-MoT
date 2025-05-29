def rotate_matrix_90_clockwise(matrix):
    # Number of rows and columns
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    
    # Create a new matrix for the rotated result
    rotated_matrix = [[0] * num_rows for _ in range(num_cols)]
    
    # Fill the rotated matrix
    for r in range(num_rows):
        for c in range(num_cols):
            rotated_matrix[c][num_rows - 1 - r] = matrix[r][c]
    
    return rotated_matrix

# Original matrix
matrix = [
    [8, 3, 2],
    [6, 3, 1]
]

# Rotate the matrix
rotated = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated)