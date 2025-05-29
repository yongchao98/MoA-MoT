# Initial matrix
matrix = [
    [8, 32],
    [0, 16]
]

# Move 'd' (down)
matrix[1][0] = matrix[0][0]  # Move 8 down
matrix[0][0] = 0  # Clear the original position

# Move 'r' (right)
matrix[1][1] = matrix[1][0]  # Move 8 right
matrix[1][0] = 0  # Clear the original position

# Print the final matrix
print(matrix)