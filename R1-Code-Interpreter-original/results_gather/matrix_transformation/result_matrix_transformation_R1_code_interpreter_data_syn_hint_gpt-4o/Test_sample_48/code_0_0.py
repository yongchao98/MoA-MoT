# Original matrix
matrix = [
    ['F', 'M', 'R'],
    ['M', 'Q', 'W'],
    ['H', 'S', 'C']
]

# Number of rows (or columns, since it's a square matrix)
n = len(matrix)

# Initialize the rotated matrix with empty values
rotated_matrix = [['' for _ in range(n)] for _ in range(n)]

# Fill the rotated matrix according to the transformation rule
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix to a single list of values
flattened_rotated_matrix = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(f"<<<{','.join(flattened_rotated_matrix)}>>>")