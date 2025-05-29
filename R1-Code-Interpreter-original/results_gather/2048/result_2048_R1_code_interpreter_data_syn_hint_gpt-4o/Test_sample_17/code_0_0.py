# Initial grid
grid = [
    [4, 0],
    [0, 16]
]

# Move sequence
moves = 'rr'

# Function to move right
def move_right(grid):
    for i in range(len(grid)):
        if grid[i][1] == 0 and grid[i][0] != 0:
            grid[i][1] = grid[i][0]
            grid[i][0] = 0
    return grid

# Apply moves
for move in moves:
    if move == 'r':
        grid = move_right(grid)

# Print the final grid
print(grid)