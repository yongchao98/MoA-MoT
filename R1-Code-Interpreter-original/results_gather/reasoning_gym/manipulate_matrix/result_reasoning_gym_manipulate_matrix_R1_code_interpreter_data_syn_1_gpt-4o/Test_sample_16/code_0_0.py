# Original matrix
matrix = [
    [2, 8, 8, 4, 8, 3, 6],
    [2, 2, 1, 7, 3, 3, 8],
    [8, 8, 2, 7, 3, 0, 1],
    [6, 1, 2, 7, 2, 2, 0],
    [8, 4, 5, 5, 1, 3, 4],
    [4, 0, 9, 4, 2, 5, 6],
    [2, 3, 6, 9, 4, 0, 5],
    [1, 4, 9, 4, 9, 7, 1],
    [1, 3, 2, 5, 7, 8, 5],
    [6, 1, 0, 5, 6, 6, 5]
]

# Cropping the matrix to rows 7-10 and columns 4-6 (1-indexed)
# Which corresponds to rows 6-9 and columns 3-5 (0-indexed)
cropped_matrix = [row[3:6] for row in matrix[6:10]]

# Print the final matrix
for row in cropped_matrix:
    print(' '.join(map(str, row)))