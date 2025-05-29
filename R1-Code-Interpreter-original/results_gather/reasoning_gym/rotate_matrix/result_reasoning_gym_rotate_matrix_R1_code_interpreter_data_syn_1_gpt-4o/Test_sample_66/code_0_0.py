# Original matrix
matrix = [
    [9, 0, 5, 6, 0, 7, 4],
    [8, 1, 7, 3, 6, 0, 5],
    [0, 0, 6, 4, 5, 0, 7],
    [3, 8, 2, 6, 2, 9, 8],
    [3, 4, 9, 8, 2, 3, 8],
    [3, 6, 8, 1, 8, 2, 6],
    [2, 5, 7, 2, 4, 8, 1]
]

# Rotate the matrix by 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)