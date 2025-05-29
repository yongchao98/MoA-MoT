# Define the matrix
matrix = [
    [6, 3, 9, 2, 3, 3, 2, 1],
    [4, 6, 8, 8, 1, 9, 7, 3],
    [6, 9, 9, 3, 5, 6, 3, 2],
    [7, 2, 9, 5, 1, 8, 1, 7],
    [5, 6, 0, 4, 8, 0, 6, 3],
    [5, 8, 1, 7, 8, 3, 9, 7],
    [0, 5, 2, 3, 2, 7, 9, 2],
    [8, 2, 0, 1, 0, 2, 7, 6]
]

# Transpose the matrix
n = len(matrix)
for i in range(n):
    for j in range(i, n):
        matrix[i][j], matrix[j][i] = matrix[j][i], matrix[i][j]

# Reverse each row
for i in range(n):
    matrix[i].reverse()

# Print the rotated matrix
print(matrix)