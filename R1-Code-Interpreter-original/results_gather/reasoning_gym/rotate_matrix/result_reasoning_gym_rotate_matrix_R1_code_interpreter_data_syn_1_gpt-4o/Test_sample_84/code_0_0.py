def rotate_90_clockwise(matrix):
    n = len(matrix)
    rotated_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rotated_matrix[j][n - 1 - i] = matrix[i][j]
    return rotated_matrix

# Input matrix
matrix = [
    [4, 0, 1, 2, 7, 7, 3, 4, 6, 9],
    [6, 9, 1, 7, 1, 5, 8, 5, 7, 0],
    [8, 8, 1, 6, 7, 2, 8, 5, 2, 7],
    [8, 1, 8, 7, 9, 3, 9, 6, 8, 4],
    [5, 6, 0, 6, 0, 6, 3, 5, 9, 1],
    [0, 6, 5, 4, 3, 5, 7, 6, 2, 4],
    [0, 6, 6, 1, 0, 4, 0, 2, 5, 4],
    [0, 7, 0, 2, 4, 8, 7, 0, 9, 9],
    [2, 7, 3, 8, 7, 6, 8, 7, 6, 8],
    [4, 5, 6, 9, 0, 1, 3, 4, 5, 8]
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
for row in rotated_matrix:
    print(row)