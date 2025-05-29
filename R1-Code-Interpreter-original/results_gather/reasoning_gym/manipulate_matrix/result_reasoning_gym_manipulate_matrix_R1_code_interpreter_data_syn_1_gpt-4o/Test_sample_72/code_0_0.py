# Initial matrix
matrix = [
    [3, 6, 5, 0, 2],
    [0, 2, 1, 2, 1]
]

# Step 2: Set all elements divisible by 6 to zero
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] % 6 == 0:
            matrix[i][j] = 0

# Step 3: Mirror the matrix along the diagonal (transpose for non-square matrix)
transposed_matrix = [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

# Print the final matrix
for row in transposed_matrix:
    print(' '.join(map(str, row)))