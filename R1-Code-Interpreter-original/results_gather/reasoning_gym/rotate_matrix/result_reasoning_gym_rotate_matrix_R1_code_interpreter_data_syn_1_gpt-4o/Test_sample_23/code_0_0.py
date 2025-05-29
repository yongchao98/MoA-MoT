def rotate_matrix_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix to store the rotated result
    rotated_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rotated_matrix[j][n - 1 - i] = matrix[i][j]
    return rotated_matrix

# Original matrix
matrix = [
    [3, 2, 4],
    [1, 2, 4],
    [1, 5, 7]
]

# Rotate the matrix by 90 degrees clockwise four times (360 degrees)
for _ in range(4):
    matrix = rotate_matrix_90_clockwise(matrix)

# Print the final matrix
print(matrix)