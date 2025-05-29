# Initial matrix
matrix = [
    [4, 3],
    [3, 2],
    [1, 8]
]

# Rotate the matrix 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the final matrix
for row in rotated_matrix:
    print(' '.join(map(str, row)))