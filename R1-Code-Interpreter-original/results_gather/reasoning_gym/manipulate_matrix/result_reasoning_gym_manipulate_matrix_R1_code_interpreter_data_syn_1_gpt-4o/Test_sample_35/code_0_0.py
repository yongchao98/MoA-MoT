# Initial matrix
matrix = [
    [6, 7, 6],
    [4, 9, 6],
    [7, 9, 4]
]

# Step 2: Set all elements divisible by 3 to zero
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] % 3 == 0:
            matrix[i][j] = 0

# Step 3: Vertically mirror the matrix
for i in range(len(matrix)):
    matrix[i][0], matrix[i][2] = matrix[i][2], matrix[i][0]

# Step 4: Horizontally mirror the matrix
matrix[0], matrix[2] = matrix[2], matrix[0]

# Step 5: Horizontally mirror the matrix again
matrix[0], matrix[2] = matrix[2], matrix[0]

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))