# Initial matrix
matrix = [
    [0, 16],
    [4, 2]
]

# Move 'd' (down)
# The number 2 is already at the bottom row, so no change is needed for the 'd' move.

# Move 'l' (left)
# Swap the numbers in the second row
matrix[1][0], matrix[1][1] = matrix[1][1], matrix[1][0]

print(matrix)