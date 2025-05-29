def rotate_90_clockwise(matrix):
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
    [3, 7, 0],
    [0, 1, 5]
]

# Rotate the matrix by 90 degrees clockwise twice
rotated_once = rotate_90_clockwise(matrix)
rotated_twice = rotate_90_clockwise(rotated_once)

# Print the result after two 90-degree rotations (180 degrees)
print(rotated_twice)