# Initial matrix
matrix = [
    [6, 3, 7, 4, 4, 8],
    [1, 7, 2, 7, 8, 2],
    [1, 5, 2, 6, 7, 8]
]

# Perform the mapping of 6 to 2
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] == 6:
            matrix[i][j] = 2

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))