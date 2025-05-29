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
    ["H", "C"],
    ["R", "U"]
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Format the output as a list of values and symbols
output = []
for row in rotated_matrix:
    output.extend(row)

print(output)