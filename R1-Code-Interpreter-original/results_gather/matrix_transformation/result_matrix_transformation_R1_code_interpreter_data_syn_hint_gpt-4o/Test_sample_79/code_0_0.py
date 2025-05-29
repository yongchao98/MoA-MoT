# Original matrix
matrix = [
    ['X', '5', 'X', 'H'],
    ['1', '9', 'W', 'K'],
    ['4', 'Q', 'K', 'L']
]

# Dimensions of the original matrix
m = len(matrix)
n = len(matrix[0])

# Initialize the rotated matrix with dimensions n x m
rotated_matrix = [['' for _ in range(m)] for _ in range(n)]

# Fill the rotated matrix
for i in range(m):
    for j in range(n):
        rotated_matrix[j][m-1-i] = matrix[i][j]

# Flatten the rotated matrix to a single list
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(flattened_result)