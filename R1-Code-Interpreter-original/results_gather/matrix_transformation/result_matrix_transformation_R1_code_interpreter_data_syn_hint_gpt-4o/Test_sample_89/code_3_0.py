def rotate_90_clockwise(matrix):
    # Determine the number of rows and columns
    num_rows = len(matrix)
    num_cols = max(len(row) for row in matrix)
    
    # Create a new matrix for the rotated version
    rotated_matrix = []
    
    # Fill the rotated matrix
    for j in range(num_cols):
        new_row = []
        for i in range(num_rows - 1, -1, -1):
            if j < len(matrix[i]):
                new_row.append(matrix[i][j])
        rotated_matrix.append(new_row)
    
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

# Format the output as a list of values and symbols
formatted_output = ','.join([','.join(row) for row in rotated_matrix])

# Print the formatted output
print(f"<<<{formatted_output}>>>")