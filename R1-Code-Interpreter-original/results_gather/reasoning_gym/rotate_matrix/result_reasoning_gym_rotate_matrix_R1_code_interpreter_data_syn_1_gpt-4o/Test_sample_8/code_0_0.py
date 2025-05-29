def rotate_matrix_90_clockwise(matrix):
    n = len(matrix)
    rotated_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rotated_matrix[j][n - 1 - i] = matrix[i][j]
    return rotated_matrix

# Input matrix
matrix = [
    [9, 8, 0, 3, 4, 3, 2, 5, 4],
    [7, 0, 8, 2, 0, 3, 9, 7, 0],
    [2, 9, 5, 1, 4, 5, 2, 8, 8],
    [1, 6, 0, 5, 5, 3, 2, 4, 2],
    [8, 4, 4, 7, 2, 9, 3, 7, 6],
    [8, 1, 2, 9, 3, 9, 9, 0, 1],
    [7, 9, 4, 1, 3, 3, 7, 7, 1],
    [7, 5, 9, 7, 4, 8, 9, 3, 1],
    [0, 2, 9, 6, 2, 7, 1, 6, 3]
]

# Rotate the matrix by 90 degrees clockwise
rotated_matrix = rotate_matrix_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)