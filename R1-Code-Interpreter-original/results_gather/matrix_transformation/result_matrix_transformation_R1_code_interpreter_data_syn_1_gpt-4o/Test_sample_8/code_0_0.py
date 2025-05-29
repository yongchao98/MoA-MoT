# Original matrix
matrix = [
    ['U', 'M'],
    ['P', 'P']
]

# Dimensions of the matrix
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix with swapped dimensions
rotated_matrix = [['' for _ in range(rows)] for _ in range(cols)]

# Fill the rotated matrix
for r in range(rows):
    for c in range(cols):
        rotated_matrix[c][rows - 1 - r] = matrix[r][c]

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)