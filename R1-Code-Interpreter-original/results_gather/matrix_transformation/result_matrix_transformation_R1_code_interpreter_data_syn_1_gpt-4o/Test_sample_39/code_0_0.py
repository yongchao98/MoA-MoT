def rotate_matrix_90_clockwise(matrix):
    # Number of rows and columns
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    
    # Create a new matrix for the rotated result
    rotated_matrix = [[None] * num_rows for _ in range(num_cols)]
    
    # Fill the rotated matrix
    for r in range(num_rows):
        for c in range(num_cols):
            rotated_matrix[c][num_rows - 1 - r] = matrix[r][c]
    
    return rotated_matrix

# Given matrix
matrix = [
    ["5", "9", "4"],
    ["8", "2", "1"]
]

# Rotate the matrix
rotated = rotate_matrix_90_clockwise(matrix)

# Format the output as a list of values and symbols
output = []
for row in rotated:
    output.extend(row)

print(output)