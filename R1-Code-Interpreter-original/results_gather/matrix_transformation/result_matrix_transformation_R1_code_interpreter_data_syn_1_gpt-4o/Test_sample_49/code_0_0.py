# Original matrix
matrix = [
    ['Q', 'P', 'K'],
    ['B', 'Z', 'I'],
    ['P', 'H', 'Y']
]

# Number of rows and columns
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix
rotated_matrix = []

# Perform the rotation
for col in range(cols):
    new_row = []
    for row in range(rows-1, -1, -1):
        new_row.append(matrix[row][col])
    rotated_matrix.append(new_row)

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)