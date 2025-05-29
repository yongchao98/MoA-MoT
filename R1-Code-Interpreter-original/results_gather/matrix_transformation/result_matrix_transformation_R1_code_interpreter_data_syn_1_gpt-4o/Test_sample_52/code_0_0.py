# Original matrix
matrix = [
    ['B', 'S', 'N'],
    ['S', 'O', 'B'],
    ['I', 'P', 'I']
]

# Number of rows (or columns, since it's a square matrix)
n = len(matrix)

# Initialize the rotated matrix with empty strings
rotated_matrix = [['' for _ in range(n)] for _ in range(n)]

# Fill the rotated matrix according to the transformation rule
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix to a single list for the output
flattened_rotated_matrix = [item for row in rotated_matrix for item in row]

# Print the result
print(flattened_rotated_matrix)