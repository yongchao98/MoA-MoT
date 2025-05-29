# Define the input matrix
matrix = [
    [3, 3, 0, 3, 7],
    [6, 5, 4, 4, 7],
    [9, 6, 4, 9, 1],
    [0, 7, 3, 7, 1],
    [7, 8, 6, 5, 0]
]

# Get the size of the matrix
n = len(matrix)

# Create a new matrix for the rotated result
rotated_matrix = [[0] * n for _ in range(n)]

# Rotate the matrix by 90 degrees clockwise
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n - 1 - i] = matrix[i][j]

# Print the rotated matrix
print(rotated_matrix)