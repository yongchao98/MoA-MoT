# Original matrix
matrix = [
    ['C', 'Z'],
    ['I', 'U']
]

# Number of rows and columns
rows = len(matrix)
cols = len(matrix[0])

# Create a new matrix for the rotated result
rotated_matrix = [['' for _ in range(rows)] for _ in range(cols)]

# Perform the rotation
for r in range(rows):
    for c in range(cols):
        rotated_matrix[c][rows - 1 - r] = matrix[r][c]

# Flatten the rotated matrix to a list of strings
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)