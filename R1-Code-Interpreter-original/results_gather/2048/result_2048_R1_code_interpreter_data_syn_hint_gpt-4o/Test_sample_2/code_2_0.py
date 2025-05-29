# Initial grid
grid = [
    [16, 32],
    [8, 16]
]

# Move right ('r')
grid[0][1], grid[0][0] = grid[0][0], 0  # Move 16 to the right
# No combination occurs because 16 and 32 are different

# Move down ('d')
grid[1][1], grid[0][1] = grid[0][1], 0  # Move 16 down
# Check for combination
if grid[1][1] == grid[1][0]:
    grid[1][1] *= 2
    grid[1][0] = 0

print(grid)