def rotate_matrix_180(matrix):
    n = len(matrix)
    # Create a new matrix with reversed rows and columns
    rotated = [[matrix[n-1-i][n-1-j] for j in range(n)] for i in range(n)]
    return rotated

# Input matrix
matrix = [
    [1, 5, 1, 4, 2, 3, 5, 1, 7, 7],
    [8, 4, 7, 2, 5, 4, 8, 9, 0, 8],
    [0, 2, 3, 7, 2, 9, 6, 3, 4, 4],
    [6, 0, 8, 3, 2, 9, 0, 3, 5, 0],
    [0, 0, 1, 6, 9, 7, 6, 2, 2, 9],
    [2, 6, 7, 7, 6, 9, 1, 3, 2, 0],
    [7, 7, 2, 7, 2, 3, 4, 0, 4, 8],
    [4, 8, 4, 5, 2, 5, 1, 3, 9, 2],
    [8, 5, 6, 3, 5, 1, 1, 9, 7, 3],
    [2, 9, 4, 4, 8, 2, 0, 8, 9, 7]
]

# Rotate the matrix
result = rotate_matrix_180(matrix)

# Print the result in the required format
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())