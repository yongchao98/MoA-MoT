def rotate_90_clockwise(matrix):
    # Determine the number of rows and columns
    num_rows = len(matrix)
    num_cols = max(len(row) for row in matrix)
    
    # Create a new matrix for the rotated version
    rotated_matrix = [['' for _ in range(num_rows)] for _ in range(num_cols)]
    
    # Fill the rotated matrix
    for i in range(num_rows):
        for j in range(len(matrix[i])):
            rotated_matrix[j][num_rows - 1 - i] = matrix[i][j]
    
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
formatted_output = ','.join([','.join(filter(None, row)) for row in rotated_matrix])

# Print the formatted output
print(f"<<<{formatted_output}>>>")