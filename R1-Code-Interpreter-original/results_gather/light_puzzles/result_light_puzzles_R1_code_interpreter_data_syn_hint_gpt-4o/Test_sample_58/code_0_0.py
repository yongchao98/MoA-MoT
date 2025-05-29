# Initial grid setup
grid = [
    [1, 0],
    [0, 0]
]

# Button A controls these positions (1-indexed in the problem statement)
controlled_positions = [(0, 1), (1, 1)]  # Converted to 0-indexed

# Toggle the lights at the controlled positions
for x, y in controlled_positions:
    grid[x][y] = 1 - grid[x][y]  # Toggle operation

# Flatten the grid to a single list for the final output
result = [grid[i][j] for i in range(2) for j in range(2)]
print(result)