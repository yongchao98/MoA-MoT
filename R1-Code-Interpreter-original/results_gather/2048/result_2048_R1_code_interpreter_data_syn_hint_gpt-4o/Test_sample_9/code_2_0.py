# Initial grid
grid = [
    [8, 0],
    [4, 0]
]

# Move sequence: 'ur'
# Move up
grid[0][0] = grid[1][0]  # Move 4 from (1, 0) to (0, 0)
grid[1][0] = 0           # Clear the original position

# Move right
grid[0][1] = grid[0][0]  # Move 4 from (0, 0) to (0, 1)
grid[0][0] = 0           # Clear the original position

# Print the final grid
print(grid)