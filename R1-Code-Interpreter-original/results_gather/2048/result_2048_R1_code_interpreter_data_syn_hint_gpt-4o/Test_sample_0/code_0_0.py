# Initial matrix
matrix = [
    [0, 16],
    [4, 2]
]

# Move 'd' (down)
matrix[1][1] = matrix[0][1]  # Move 2 down
matrix[0][1] = 0  # Clear the original position

# Move 'l' (left)
if matrix[1][0] == matrix[1][1]:  # Check if they can combine
    matrix[1][0] *= 2  # Combine to form 8
    matrix[1][1] = 0  # Clear the original position
else:
    matrix[1][0], matrix[1][1] = matrix[1][1], matrix[1][0]  # Swap positions

print(matrix)