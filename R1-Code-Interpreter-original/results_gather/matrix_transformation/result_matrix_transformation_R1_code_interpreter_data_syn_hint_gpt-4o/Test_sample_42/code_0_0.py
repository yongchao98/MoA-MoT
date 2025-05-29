# Original matrix
matrix = [
    ['S', 'L', 'D'],
    ['R', 'I', 'S'],
    ['K', 'F', 'L']
]

# Number of rows (or columns, since it's a square matrix)
n = len(matrix)

# Create a new matrix for the rotated result
rotated_matrix = [['' for _ in range(n)] for _ in range(n)]

# Apply the rotation transformation
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix to a list of values
flattened_result = [item for row in rotated_matrix for item in row]

# Print the result in the required format
print(flattened_result)