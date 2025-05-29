# Original matrix
matrix = [
    [1, 9, 1, 6, 7],
    [5, 1, 7, 1, 2],
    [3, 9, 0, 0, 3],
    [1, 2, 9, 8, 4],
    [2, 6, 9, 6, 3],
    [9, 9, 3, 1, 5],
    [3, 9, 3, 9, 8],
    [7, 4, 4, 3, 8],
    [7, 8, 3, 8, 4]
]

# Dimensions of the matrix
rows = len(matrix)
cols = len(matrix[0])

# Create a new matrix for the mirrored result
mirrored_matrix = [[0] * rows for _ in range(cols)]

# Perform the mirroring along the counterdiagonal
for i in range(rows):
    for j in range(cols):
        mirrored_matrix[cols - j - 1][rows - i - 1] = matrix[i][j]

# Print the final mirrored matrix
for row in mirrored_matrix:
    print(' '.join(map(str, row)))