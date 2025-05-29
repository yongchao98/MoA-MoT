# Original matrix
matrix = [
    ['X', '4', '2', 'N'],
    ['H', '1', '7', '2'],
    ['M', '8', 'K', 'I']
]

# Number of rows and columns in the original matrix
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix with dimensions cols x rows
rotated_matrix = [['' for _ in range(rows)] for _ in range(cols)]

# Fill the rotated matrix
for r in range(rows):
    for c in range(cols):
        rotated_matrix[c][rows - 1 - r] = matrix[r][c]

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)