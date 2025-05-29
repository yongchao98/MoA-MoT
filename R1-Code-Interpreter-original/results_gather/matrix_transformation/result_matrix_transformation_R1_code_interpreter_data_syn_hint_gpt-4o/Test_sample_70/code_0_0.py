# Original matrix
matrix = [
    ['O', 'U', 'Y', '1'],
    ['T', '9', 'Q', 'W'],
    ['Y', '8', 'Z', 'D']
]

# Dimensions of the matrix
m = len(matrix)
n = len(matrix[0])

# Initialize the rotated matrix with dimensions n x m
rotated_matrix = [['' for _ in range(m)] for _ in range(n)]

# Fill the rotated matrix
for i in range(m):
    for j in range(n):
        rotated_matrix[j][m-1-i] = matrix[i][j]

# Flatten the rotated matrix to a single list
flattened_rotated_matrix = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_rotated_matrix)