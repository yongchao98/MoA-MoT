# Initial grid
grid = [
    [1, 1],
    [1, 0]
]

# Positions controlled by Button A
positions = [(1, 1), (0, 1)]

# Toggle the lights at the specified positions
for x, y in positions:
    grid[x][y] = 1 - grid[x][y]

# Flatten the grid to a single list for output
result = [grid[i][j] for i in range(2) for j in range(2)]

# Print the result
print(result)