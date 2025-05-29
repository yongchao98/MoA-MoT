# Original matrix
matrix = [
    [8, 7, 9],
    [2, 7, 4],
    [9, 1, 7],
    [6, 8, 8],
    [3, 7, 1],
    [0, 4, 6]
]

# Dimensions of the matrix
m = len(matrix)
n = len(matrix[0])

# Create a new matrix to store the result
mirrored_matrix = [[0] * m for _ in range(n)]

# Perform the mirroring along the counterdiagonal
for i in range(m):
    for j in range(n):
        mirrored_matrix[n - j - 1][m - i - 1] = matrix[i][j]

# Print the final mirrored matrix
for row in mirrored_matrix:
    print(' '.join(map(str, row)))