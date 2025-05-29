def rotate_matrix_90(matrix):
    n = len(matrix)
    # Create a new matrix with rotated values
    rotated = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rotated[j][n-1-i] = matrix[i][j]
    return rotated

# Initial matrix
matrix = [[8, 5],
          [0, 1]]

# Simulate 720 degrees rotation (8 times 90 degrees)
result = matrix
for _ in range(8):  # 8 times because 720/90 = 8
    result = rotate_matrix_90(result)

# Print the result in the required format
print(f"{result[0][0]} {result[0][1]}\n{result[1][0]} {result[1][1]}")