# Define the matrix
matrix = [
    [9, 9, 4, 7, 6, 2, 5, 3, 8],
    [7, 0, 6, 9, 3, 7, 7, 4, 5],
    [3, 3, 9, 2, 7, 8, 5, 6, 9],
    [6, 4, 0, 1, 1, 3, 4, 6, 1],
    [9, 0, 1, 3, 4, 8, 7, 0, 4],
    [6, 8, 2, 4, 8, 6, 2, 8, 9],
    [1, 3, 9, 4, 3, 9, 9, 2, 6],
    [6, 3, 9, 2, 1, 5, 9, 1, 9],
    [1, 8, 5, 1, 7, 1, 7, 8, 8]
]

# Rotate the matrix by 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)