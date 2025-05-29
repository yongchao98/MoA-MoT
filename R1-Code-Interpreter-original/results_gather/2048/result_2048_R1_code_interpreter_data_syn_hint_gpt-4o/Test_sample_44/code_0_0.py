# Initial matrix
matrix = [
    [2, 4, 32, 0],
    [0, 0, 16, 0],
    [0, 0, 0, 0],
    [2, 8, 0, 0]
]

# Initial position of the number
position = (1, 3)  # (row, column) in 0-indexed

# Move sequence
moves = 'lrud'

# Function to perform a move
def move_number(matrix, position, direction):
    row, col = position
    if direction == 'l' and col > 0:
        col -= 1
    elif direction == 'r' and col < 3:
        col += 1
    elif direction == 'u' and row > 0:
        row -= 1
    elif direction == 'd' and row < 3:
        row += 1
    return (row, col)

# Simulate the moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)