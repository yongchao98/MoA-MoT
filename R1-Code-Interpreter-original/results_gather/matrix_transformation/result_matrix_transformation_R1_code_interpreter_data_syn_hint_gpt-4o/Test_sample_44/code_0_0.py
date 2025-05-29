# Original matrix
matrix = [
    ['T', 'K', 'M'],
    ['C', 'M', 'A'],
    ['X', 'U', 'I']
]

# Number of rows (or columns, since it's a square matrix)
n = len(matrix)

# Create a new matrix for the rotated version
rotated_matrix = [['' for _ in range(n)] for _ in range(n)]

# Apply the transformation
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix to a single list
flattened_rotated_matrix = [item for row in rotated_matrix for item in row]

# Print the result in the required format
print(f"<<<{','.join(flattened_rotated_matrix)}>>>")