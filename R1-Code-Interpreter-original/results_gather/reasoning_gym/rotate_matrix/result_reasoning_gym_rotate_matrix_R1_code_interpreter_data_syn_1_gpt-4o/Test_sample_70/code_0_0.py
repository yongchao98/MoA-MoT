def rotate_matrix_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix to store the rotated result
    rotated_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rotated_matrix[j][n - 1 - i] = matrix[i][j]
    return rotated_matrix

def rotate_matrix_720_clockwise(matrix):
    # Rotate the matrix 4 times by 90 degrees to achieve 360 degrees twice (720 degrees)
    for _ in range(4):
        matrix = rotate_matrix_90_clockwise(matrix)
    return matrix

# Input matrix
matrix = [
    [9, 5, 2, 8],
    [9, 2, 6, 4],
    [9, 6, 9, 9],
    [0, 8, 8, 9]
]

# Rotate the matrix by 720 degrees
result_matrix = rotate_matrix_720_clockwise(matrix)

# Print the result
print(result_matrix)