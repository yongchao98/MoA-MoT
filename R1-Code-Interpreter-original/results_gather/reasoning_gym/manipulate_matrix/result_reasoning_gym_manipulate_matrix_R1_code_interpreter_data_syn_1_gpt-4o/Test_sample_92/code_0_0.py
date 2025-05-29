# Initial matrix
matrix = [
    [0, 2, 0, 6, 4, 5, 0, 5, 2, 6],
    [4, 6, 0, 8, 7, 9, 6, 0, 3, 3],
    [1, 6, 5, 3, 2, 6, 3, 5, 8, 4],
    [6, 2, 1, 7, 8, 9, 1, 6, 9, 8],
    [7, 1, 9, 7, 1, 2, 8, 6, 3, 2],
    [7, 8, 0, 7, 6, 5, 9, 7, 1, 2]
]

# Step 2: Map each occurrence of 0 to 9
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] == 0:
            matrix[i][j] = 9

# Step 3: Vertically mirror the matrix
for i in range(len(matrix)):
    matrix[i] = matrix[i][::-1]

# Print the final matrix
for row in matrix:
    print(" ".join(map(str, row)))