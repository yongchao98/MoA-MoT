def rotate_90_clockwise(matrix):
    n = len(matrix)
    # Create a new matrix to store the rotated result
    rotated_matrix = [[0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(n):
            rotated_matrix[j][n - 1 - i] = matrix[i][j]
    
    return rotated_matrix

# Input matrix
matrix = [
    [9, 3, 6, 7, 4, 8, 1, 8],
    [6, 8, 5, 0, 1, 8, 1, 2],
    [7, 0, 8, 6, 9, 3, 4, 4],
    [5, 5, 8, 1, 6, 8, 7, 4],
    [2, 3, 4, 6, 7, 3, 2, 7],
    [3, 1, 7, 5, 5, 3, 2, 4],
    [2, 3, 4, 2, 0, 9, 4, 2],
    [0, 7, 6, 3, 3, 7, 8, 7]
]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)