def rotate_90_clockwise(matrix):
    # Determine the maximum number of columns in any row
    max_columns = max(len(row) for row in matrix)
    
    # Create a new matrix for the rotated version
    rotated_matrix = []
    
    # Fill the rotated matrix
    for j in range(max_columns):
        new_row = []
        for i in range(len(matrix) - 1, -1, -1):
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