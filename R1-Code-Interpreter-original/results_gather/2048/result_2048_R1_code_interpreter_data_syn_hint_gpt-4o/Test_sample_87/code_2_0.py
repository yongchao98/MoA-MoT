# Initial grid
grid = [
    [0, 0, 0, 0, 0],
    [0, 0, 2, 8, 0],
    [0, 0, 0, 8, 0],
    [32, 0, 0, 0, 0],
    [0, 8, 0, 0, 0]
]

# Initial position of the number to move
position = (1, 3)  # Adjusted to 0-based index

# Sequence of moves
moves = 'rlulur'

# Function to perform a move
def move_and_combine(grid, position, direction):
    x, y = position
    if direction == 'r':
        if y < 4 and grid[x][y+1] == 0:
            grid[x][y+1] = grid[x][y]
            grid[x][y] = 0
            y += 1
        elif y < 4 and grid[x][y+1] == grid[x][y]:
            grid[x][y+1] *= 2
            grid[x][y] = 0
            y += 1
    elif direction == 'l':
        if y > 0 and grid[x][y-1] == 0:
            grid[x][y-1] = grid[x][y]
            grid[x][y] = 0
            y -= 1
        elif y > 0 and grid[x][y-1] == grid[x][y]:
            grid[x][y-1] *= 2
            grid[x][y] = 0
            y -= 1
    elif direction == 'u':
        if x > 0 and grid[x-1][y] == 0:
            grid[x-1][y] = grid[x][y]
            grid[x][y] = 0
            x -= 1
        elif x > 0 and grid[x-1][y] == grid[x][y]:
            grid[x-1][y] *= 2
            grid[x][y] = 0
            x -= 1
    elif direction == 'd':
        if x < 4 and grid[x+1][y] == 0:
            grid[x+1][y] = grid[x][y]
            grid[x][y] = 0
            x += 1
        elif x < 4 and grid[x+1][y] == grid[x][y]:
            grid[x+1][y] *= 2
            grid[x][y] = 0
            x += 1
    return (x, y)

# Process each move
for move in moves:
    position = move_and_combine(grid, position, move)

# Print the final grid
print(grid)