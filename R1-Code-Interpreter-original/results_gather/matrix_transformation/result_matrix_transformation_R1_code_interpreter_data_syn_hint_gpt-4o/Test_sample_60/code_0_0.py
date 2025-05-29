def rotate_90_clockwise(matrix):
    # Number of rows and columns
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    
    # Create a new matrix for the rotated result
    rotated_matrix = []
    
    # Rotate the matrix
    for col in range(num_cols):
        new_row = []
        for row in range(num_rows - 1, -1, -1):
            new_row.append(matrix[row][col])
        rotated_matrix.append(new_row)
    
    return rotated_matrix

# Original matrix
matrix = [
    ['A', 'J', 'X', 'X'],
    ['9', 'S', '4', '7'],
    ['R', 'X', 'V', 'T']
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)