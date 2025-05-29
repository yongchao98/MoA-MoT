# Original matrix
matrix = [
    ['9', 'H', 'B', 'X'],
    ['S', 'X', 'T', '6'],
    ['U', 'V', 'N', '6']
]

# Number of rows and columns
n = len(matrix)
m = len(matrix[0])

# Initialize the rotated matrix
rotated_matrix = [['' for _ in range(n)] for _ in range(m)]

# Perform the rotation
for i in range(n):
    for j in range(m):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)