# Original matrix
matrix = [
    [4, 0],
    [7, 8]
]

# Function to rotate the matrix 90 degrees clockwise
def rotate_90_clockwise(mat):
    n = len(mat)
    # Create a new matrix to store the rotated result
    rotated = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rotated[j][n - 1 - i] = mat[i][j]
    return rotated

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)