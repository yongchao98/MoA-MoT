# Initial grid
grid = [
    [0, 0, 0, 16, 16],
    [0, 0, 0, 2, 2],
    [0, 0, 0, 0, 0],
    [0, 0, 16, 4, 0],
    [0, 0, 0, 0, 0]
]

# Initial position
x, y = 2, 3

# Move sequence
moves = 'rlrdr'

# Function to perform a move
def move(x, y, direction):
    if direction == 'u':
        return max(x - 1, 0), y
    elif direction == 'd':
        return min(x + 1, 4), y
    elif direction == 'l':
        return x, max(y - 1, 0)
    elif direction == 'r':
        return x, min(y + 1, 4)
    return x, y

# Simulate the moves
for move_direction in moves:
    x, y = move(x, y, move_direction)

# Since the number at (2, 3) was 0, the grid remains unchanged
print(grid)