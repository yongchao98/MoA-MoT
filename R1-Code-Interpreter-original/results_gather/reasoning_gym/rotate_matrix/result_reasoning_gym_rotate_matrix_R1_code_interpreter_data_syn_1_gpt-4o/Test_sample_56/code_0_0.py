def rotate_matrix_90_clockwise(matrix):
    n = len(matrix)
    return [[matrix[n - j - 1][i] for j in range(n)] for i in range(n)]

def rotate_matrix_720_clockwise(matrix):
    # 720 degrees is equivalent to two full rotations (360 degrees each)
    # which means the matrix will be the same as the original
    return matrix

matrix = [
    [3, 5, 2, 0, 8, 0, 5, 7],
    [2, 3, 1, 2, 9, 2, 7, 6],
    [5, 6, 2, 4, 2, 1, 7, 5],
    [0, 2, 3, 7, 3, 3, 2, 0],
    [3, 8, 7, 8, 0, 7, 6, 8],
    [6, 7, 8, 1, 7, 2, 6, 0],
    [2, 2, 6, 4, 0, 2, 1, 8],
    [8, 1, 6, 4, 5, 3, 6, 5]
]

rotated_matrix = rotate_matrix_720_clockwise(matrix)
print(rotated_matrix)