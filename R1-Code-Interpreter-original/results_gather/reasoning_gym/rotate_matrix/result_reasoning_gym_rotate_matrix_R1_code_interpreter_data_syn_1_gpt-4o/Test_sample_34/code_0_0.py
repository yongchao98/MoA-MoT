# Define the original matrix
matrix = [
    [4, 6, 3],
    [7, 4, 4],
    [8, 1, 7]
]

# Rotate the matrix by 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)