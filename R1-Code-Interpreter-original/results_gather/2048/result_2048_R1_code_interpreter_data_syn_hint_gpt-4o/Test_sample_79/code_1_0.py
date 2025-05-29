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
    new_x, new_y = move(x, y, move_direction)
    if grid[new_x][new_y] == grid[x][y] and grid[x][y] != 0:
        # Combine numbers
        grid[new_x][new_y] *= 2
        grid[x][y] = 0
    elif grid[new_x][new_y] == 0:
        # Move number
        grid[new_x][new_y] = grid[x][y]
        grid[x][y] = 0
    x, y = new_x, new_y

print(grid)