# Original matrix
matrix = [
    ['O', 'T', '4', 'E'],
    ['O', '4', '2', 'K'],
    ['D', 'Q', 'W', 'Z']
]

# Number of rows and columns in the original matrix
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix with dimensions cols x rows
rotated_matrix = [['' for _ in range(rows)] for _ in range(cols)]

# Fill the rotated matrix using the transformation rule
for i in range(rows):
    for j in range(cols):
        rotated_matrix[j][rows - 1 - i] = matrix[i][j]

# Flatten the rotated matrix to a single list for output
flattened_rotated_matrix = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_rotated_matrix)