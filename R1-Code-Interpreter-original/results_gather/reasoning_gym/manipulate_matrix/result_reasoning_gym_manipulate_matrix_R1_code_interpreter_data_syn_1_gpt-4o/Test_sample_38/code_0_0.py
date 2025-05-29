# Original matrix
matrix = [
    [7, 7, 2, 5, 1, 2],
    [9, 3, 6, 6, 3, 6],
    [1, 8, 5, 4, 6, 7]
]

# Transpose the matrix
transposed_matrix = [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

# Print the transposed matrix
for row in transposed_matrix:
    print(' '.join(map(str, row)))