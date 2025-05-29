def rotate_270_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = [list(row) for row in zip(*matrix)]
    # Reverse the order of the rows
    rotated_matrix = transposed_matrix[::-1]
    return rotated_matrix

# Define the matrix
matrix = [
    [1, 4, 4, 6, 0],
    [2, 6, 3, 7, 6],
    [8, 0, 3, 2, 5],
    [6, 5, 6, 2, 2],
    [5, 5, 0, 3, 1]
]

# Rotate the matrix
rotated_matrix = rotate_270_clockwise(matrix)

# Print the rotated matrix
for row in rotated_matrix:
    print(row)