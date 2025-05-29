# Original matrix
matrix = [
    [4, 1],
    [6, 4]
]

# Size of the matrix
n = len(matrix)

# Create a new matrix to store the rotated result
rotated_matrix = [[0] * n for _ in range(n)]

# Perform the rotation
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n - 1 - i] = matrix[i][j]

# Print the rotated matrix
print(rotated_matrix)