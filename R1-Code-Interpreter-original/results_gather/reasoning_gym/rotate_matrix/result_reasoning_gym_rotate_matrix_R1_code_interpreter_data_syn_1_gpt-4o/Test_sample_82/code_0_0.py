def rotate_90_clockwise(matrix):
    n = len(matrix)
    return [[matrix[n - j - 1][i] for j in range(n)] for i in range(n)]

# Original matrix
matrix = [
    [9, 1, 5, 3, 9],
    [0, 2, 8, 8, 8],
    [4, 9, 1, 8, 5],
    [5, 1, 1, 7, 6],
    [1, 9, 8, 4, 5]
]

# Rotate the matrix 4 times by 90 degrees
for _ in range(4):
    matrix = rotate_90_clockwise(matrix)

print(matrix)