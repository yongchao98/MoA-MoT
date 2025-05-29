# Define the initial matrix
matrix = [
    [2, 4, 3, 7, 0],
    [0, 1, 5, 8, 5],
    [0, 8, 7, 5, 9]
]

# Perform the mapping of 2 to 8
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] == 2:
            matrix[i][j] = 8

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))