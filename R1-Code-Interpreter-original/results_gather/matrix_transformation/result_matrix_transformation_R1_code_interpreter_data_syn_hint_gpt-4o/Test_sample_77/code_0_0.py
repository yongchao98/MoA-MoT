# Original matrix
matrix = [
    ['X', '4', '2', 'N'],
    ['H', '1', '7', '2'],
    ['M', '8', 'K', 'I']
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

# Flatten the rotated matrix to a list of values and symbols
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(f"<<<{','.join(flattened_result)}>>>")