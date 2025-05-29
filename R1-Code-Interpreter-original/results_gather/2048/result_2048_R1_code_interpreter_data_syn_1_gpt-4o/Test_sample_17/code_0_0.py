# Initial grid
grid = [
    [4, 0],
    [0, 16]
]

# Sequence of moves
moves = 'rr'

# Function to move right
def move_right(grid):
    for i in range(len(grid)):
        if grid[i][1] == 0 and grid[i][0] != 0:
            grid[i][1] = grid[i][0]
            grid[i][0] = 0
        elif grid[i][1] == grid[i][0] and grid[i][0] != 0:
            grid[i][1] *= 2
            grid[i][0] = 0
    return grid

# Apply the moves
for move in moves:
    if move == 'r':
        grid = move_right(grid)

# Print the final grid
print(grid)