# Original matrix
matrix = [
    [8, 9, 8, 8, 9, 4, 1, 9],
    [8, 7, 2, 1, 9, 5, 0, 0],
    [0, 2, 9, 2, 9, 9, 2, 6],
    [7, 2, 1, 2, 7, 1, 6, 5],
    [6, 0, 4, 2, 2, 3, 9, 2],
    [3, 9, 2, 7, 0, 6, 4, 4],
    [4, 1, 9, 8, 0, 5, 6, 2],
    [9, 9, 4, 5, 1, 1, 7, 7]
]

# Rotate the matrix by 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)