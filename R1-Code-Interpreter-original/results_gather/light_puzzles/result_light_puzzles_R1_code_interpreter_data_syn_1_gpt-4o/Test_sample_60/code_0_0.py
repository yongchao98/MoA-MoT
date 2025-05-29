# Initial grid setup
grid = [[0, 1], [0, 1]]

# Button A controls these positions (0-indexed)
positions = [(0, 0), (1, 1)]

# Toggle the lights at the specified positions
for x, y in positions:
    grid[x][y] = 1 - grid[x][y]  # Toggle operation

# Flatten the grid to a single list
result = [grid[i][j] for i in range(2) for j in range(2)]

# Print the result
print(result)