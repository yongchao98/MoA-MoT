# Original matrix
matrix = [
    [5, 6, 1],
    [5, 8, 8],
    [0, 0, 4]
]

# Rotate the matrix by 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)