# Define the initial matrix
matrix = [
    [5, 1, 4],
    [2, 3, 5],
    [1, 7, 7],
    [8, 4, 7],
    [2, 5, 4],
    [8, 9, 0],
    [8, 0, 2],
    [3, 7, 2],
    [9, 6, 3],
    [4, 4, 6]
]

# Perform the transformation: map 3 to 2
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] == 3:
            matrix[i][j] = 2

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))