# Define the input matrix
matrix = [
    [9, 9, 6, 7, 4],
    [2, 6, 8, 1, 4],
    [4, 9, 0, 4, 5],
    [3, 2, 9, 3, 2],
    [8, 2, 3, 5, 1]
]

# Determine the size of the matrix
n = len(matrix)

# Create a new matrix to store the rotated result
rotated_matrix = [[0] * n for _ in range(n)]

# Apply the rotation transformation
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n - 1 - i] = matrix[i][j]

# Print the rotated matrix
print(rotated_matrix)