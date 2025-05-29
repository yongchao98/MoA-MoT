# Define the initial matrix
matrix = [
    [8, 7],
    [5, 6],
    [0, 2],
    [2, 6],
    [2, 5],
    [7, 0]
]

# Perform the operation of mapping 7 to 0
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] == 7:
            matrix[i][j] = 0

# Print the final matrix
for row in matrix:
    print(' '.join(map(str, row)))