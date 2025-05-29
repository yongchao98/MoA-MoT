# Original matrix
matrix = [
    [8, 4, 7, 7, 6, 3, 7, 5, 9],
    [6, 3, 3, 0, 7, 8, 7, 9, 1],
    [8, 5, 8, 4, 6, 4, 8, 2, 1],
    [4, 4, 7, 9, 0, 6, 0, 5, 4]
]

# Vertically mirror the matrix
mirrored_matrix = [row[::-1] for row in matrix]

# Print the final matrix
for row in mirrored_matrix:
    print(' '.join(map(str, row)))