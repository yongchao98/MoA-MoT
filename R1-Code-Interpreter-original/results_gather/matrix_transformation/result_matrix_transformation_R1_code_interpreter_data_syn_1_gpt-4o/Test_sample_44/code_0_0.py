# Original matrix
matrix = [
    ['T', 'K', 'M'],
    ['C', 'M', 'A'],
    ['X', 'U', 'I']
]

# Dimensions of the matrix
n = len(matrix)

# Create a new matrix for the rotated result
rotated_matrix = [['' for _ in range(n)] for _ in range(n)]

# Perform the rotation
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n - 1 - i] = matrix[i][j]

# Flatten the rotated matrix to a single list
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)