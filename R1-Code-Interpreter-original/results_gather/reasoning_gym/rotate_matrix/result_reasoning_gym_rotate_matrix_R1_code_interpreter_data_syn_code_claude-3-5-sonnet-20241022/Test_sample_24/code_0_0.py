def print_matrix(matrix):
    for row in matrix:
        print(' '.join(map(str, row)))

# Input matrix
matrix = [
    [0, 7, 8, 4, 8, 6, 3, 1, 5, 9],
    [1, 5, 4, 9, 3, 0, 9, 3, 0, 9],
    [8, 2, 6, 0, 2, 6, 3, 2, 1, 9],
    [3, 8, 3, 7, 2, 6, 2, 1, 3, 8],
    [2, 3, 2, 2, 8, 7, 5, 0, 0, 9],
    [5, 0, 9, 0, 8, 8, 1, 9, 7, 6],
    [5, 7, 5, 3, 8, 7, 8, 7, 0, 3],
    [2, 5, 7, 2, 0, 3, 4, 2, 9, 4],
    [7, 0, 2, 0, 9, 4, 8, 2, 8, 1],
    [0, 9, 4, 9, 6, 3, 0, 8, 7, 8]
]

# Since we're rotating 360 degrees (4 times 90 degrees), 
# the output will be the same as input
print_matrix(matrix)