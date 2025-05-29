# Original matrix
matrix = [
    ['N', 'H', 'Y', '9'],
    ['X', 'L', 'W', 'F'],
    ['T', '1', '5', 'F']
]

# Number of rows and columns
rows = len(matrix)
cols = len(matrix[0])

# Create a new matrix for the rotated result
rotated_matrix = []

# Rotate the matrix 90 degrees clockwise
for col in range(cols):
    new_row = []
    for row in range(rows-1, -1, -1):
        new_row.append(matrix[row][col])
    rotated_matrix.append(new_row)

# Flatten the rotated matrix for the output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)