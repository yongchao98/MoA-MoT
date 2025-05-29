# Original matrix
matrix = [
    [3, 9, 6],
    [2, 9, 8],
    [2, 5, 4]
]

# Rotate the matrix by 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)