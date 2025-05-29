# Define the input matrix
matrix = [
    [7, 8, 4, 2],
    [2, 4, 2, 0],
    [6, 7, 9, 1],
    [3, 4, 8, 2]
]

# Determine the size of the matrix
n = len(matrix)

# Create a new matrix to store the rotated result
rotated_matrix = [[0] * n for _ in range(n)]

# Apply the transformation to rotate the matrix
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n - 1 - i] = matrix[i][j]

# Print the rotated matrix
print(rotated_matrix)