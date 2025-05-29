# Initial grid as a list
grid = [0, 1, 0, 1]

# Positions controlled by Button A
positions = [0, 3]

# Toggle the lights at the specified positions
for pos in positions:
    grid[pos] = 1 - grid[pos]  # Toggle between 0 and 1

# Print the final state of the grid
print(grid)