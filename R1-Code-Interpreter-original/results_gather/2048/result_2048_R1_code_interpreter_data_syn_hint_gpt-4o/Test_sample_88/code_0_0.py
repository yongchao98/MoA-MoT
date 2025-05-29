# Initial matrix
matrix = [
    [0, 4, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 32, 0, 0, 0],
    [0, 0, 0, 16, 8],
    [0, 0, 0, 16, 16]
]

# Initial position of the number to move
position = (0, 4)  # (row, column) for (1, 5) in 0-indexed

# Sequence of moves
moves = 'lrruru'

# Function to perform a move
def move(position, direction):
    row, col = position
    if direction == 'u':
        row = max(0, row - 1)
    elif direction == 'd':
        row = min(4, row + 1)
    elif direction == 'l':
        col = max(0, col - 1)
    elif direction == 'r':
        col = min(4, col + 1)
    return (row, col)

# Simulate the moves
for move_direction in moves:
    position = move(position, move_direction)

# Print the final matrix
print(matrix)