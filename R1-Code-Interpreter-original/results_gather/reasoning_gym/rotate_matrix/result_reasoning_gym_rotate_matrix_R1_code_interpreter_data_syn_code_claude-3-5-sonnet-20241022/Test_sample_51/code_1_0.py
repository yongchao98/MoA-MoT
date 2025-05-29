def rotate_180(matrix):
    n = len(matrix)
    # For 180 degrees, we reverse rows and columns
    return [[matrix[n-1-i][n-1-j] for j in range(n)] for i in range(n)]

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

# 900 degrees = 2 full rotations (720°) + 180°
# Two full rotations bring us back to original
result = rotate_180(matrix)

# Print result
for row in result:
    print(" ".join(map(str, row)))