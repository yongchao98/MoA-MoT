# Original matrix
matrix = [
    [8, 1, 4, 5, 8],
    [5, 1, 5, 0, 4],
    [6, 1, 7, 8, 9],
    [1, 3, 4, 2, 8]
]

# Step 2: Transpose the matrix
transposed_matrix = [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

# Step 3: Set elements divisible by 7 to zero
for i in range(len(transposed_matrix)):
    for j in range(len(transposed_matrix[0])):
        if transposed_matrix[i][j] % 7 == 0:
            transposed_matrix[i][j] = 0

# Print the final matrix
for row in transposed_matrix:
    print(' '.join(map(str, row)))