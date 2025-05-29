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

# Print the result in the required format
formatted_output = ','.join([','.join(row) for row in rotated_matrix])
print(f"<<<{formatted_output}>>>")