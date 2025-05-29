# Initial matrix
matrix = [
    [0, 16, 32],
    [0, 0, 32],
    [16, 4, 2]
]

# Initial position of the number
row, col = 2, 1  # 0-indexed position (3, 2)

# Move sequence
moves = 'dlr'

# Execute the moves
for move in moves:
    if move == 'd':
        if row < 2:  # Check if we can move down
            row += 1
    elif move == 'l':
        if col > 0:  # Check if we can move left
            col -= 1
    elif move == 'r':
        if col < 2:  # Check if we can move right
            col += 1

# Update the matrix with the final position of the number
matrix[2][1] = 0  # Clear the original position
matrix[row][col] = 4  # Place the number in the new position

# Print the final matrix
print(matrix)