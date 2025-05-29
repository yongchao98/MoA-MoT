def rotate_matrix_90_clockwise(matrix):
    n = len(matrix)
    return [[matrix[n - 1 - j][i] for j in range(n)] for i in range(n)]

def rotate_matrix_900_clockwise(matrix):
    # 900 degrees is equivalent to 10 * 90 degrees
    # Since 360 degrees is a full rotation, 900 degrees results in the same matrix
    return matrix

# Input matrix
matrix = [
    [7, 1, 2, 6, 0, 0, 3, 3, 2],
    [4, 1, 0, 6, 2, 5, 4, 8, 3],
    [7, 3, 7, 2, 0, 0, 0, 7, 5],
    [6, 9, 2, 8, 5, 7, 7, 9, 5],
    [7, 3, 2, 6, 8, 8, 2, 9, 2],
    [0, 5, 7, 3, 7, 9, 0, 1, 7],
    [9, 8, 6, 2, 3, 4, 4, 1, 5],
    [1, 2, 1, 9, 6, 9, 9, 1, 2],
    [4, 7, 4, 5, 6, 7, 2, 9, 5]
]

# Rotate the matrix by 900 degrees clockwise
result_matrix = rotate_matrix_900_clockwise(matrix)

# Print the result
print(result_matrix)