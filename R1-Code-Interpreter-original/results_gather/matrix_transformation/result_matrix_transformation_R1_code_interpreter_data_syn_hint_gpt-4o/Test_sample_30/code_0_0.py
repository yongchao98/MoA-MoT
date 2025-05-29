# Original matrix
matrix = [
    [6, 1, 7],
    [9, 5, 1]
]

# Number of rows and columns in the original matrix
m = len(matrix)
n = len(matrix[0])

# Initialize the rotated matrix with dimensions n x m
rotated_matrix = [[0] * m for _ in range(n)]

# Fill the rotated matrix
for i in range(m):
    for j in range(n):
        rotated_matrix[j][m - 1 - i] = matrix[i][j]

# Print the rotated matrix
print(rotated_matrix)