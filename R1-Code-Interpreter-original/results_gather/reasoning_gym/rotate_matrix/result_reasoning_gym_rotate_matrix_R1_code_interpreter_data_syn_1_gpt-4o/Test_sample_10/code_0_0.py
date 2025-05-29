# Define the original matrix
matrix = [
    [0, 4, 3],
    [3, 2, 1],
    [8, 1, 9]
]

# Rotate the matrix 180 degrees
rotated_matrix = [row[::-1] for row in matrix[::-1]]

# Print the rotated matrix
print(rotated_matrix)