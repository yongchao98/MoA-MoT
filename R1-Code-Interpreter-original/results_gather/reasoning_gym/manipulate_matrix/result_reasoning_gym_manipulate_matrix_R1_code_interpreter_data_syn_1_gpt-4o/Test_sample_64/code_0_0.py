# Original matrix
matrix = [
    [6, 3, 1, 5, 0, 0, 2],
    [6, 5, 1, 4, 5, 5, 4],
    [2, 8, 0, 9, 7, 8, 7],
    [7, 5, 3, 5, 8, 8, 5],
    [0, 2, 0, 1, 7, 8, 3],
    [5, 1, 2, 9, 1, 8, 4],
    [9, 4, 4, 3, 6, 6, 9]
]

# Vertically mirror the matrix
mirrored_matrix = [row[::-1] for row in matrix]

# Print the final matrix
for row in mirrored_matrix:
    print(' '.join(map(str, row)))