# Define the matrix
matrix = [
    [6, 8, 6, 5, 5, 8],
    [6, 5, 0, 8, 6, 3],
    [9, 0, 6, 3, 2, 8],
    [9, 0, 2, 7, 9, 2],
    [2, 3, 5, 5, 9, 2],
    [3, 4, 4, 1, 5, 1]
]

# Rotate the matrix 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)