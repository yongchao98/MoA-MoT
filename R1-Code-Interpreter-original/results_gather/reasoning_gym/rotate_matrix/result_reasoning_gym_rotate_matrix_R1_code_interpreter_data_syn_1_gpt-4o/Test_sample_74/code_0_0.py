# Define the input matrix as a list of lists
matrix = [
    [0, 8, 8, 7, 0, 4, 1, 3, 7],
    [5, 6, 7, 4, 3, 7, 4, 7, 2],
    [2, 4, 5, 7, 5, 1, 4, 5, 5],
    [2, 2, 6, 1, 0, 1, 5, 1, 2],
    [5, 7, 8, 4, 9, 4, 7, 7, 6],
    [3, 2, 2, 7, 3, 2, 3, 3, 5],
    [9, 8, 3, 2, 5, 4, 2, 9, 5],
    [8, 7, 8, 6, 3, 6, 3, 8, 5],
    [6, 6, 2, 2, 9, 0, 5, 8, 8]
]

# Get the size of the matrix
n = len(matrix)

# Create a new matrix for the rotated result
rotated_matrix = [[0] * n for _ in range(n)]

# Perform the rotation
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n - 1 - i] = matrix[i][j]

# Print the rotated matrix
print(rotated_matrix)