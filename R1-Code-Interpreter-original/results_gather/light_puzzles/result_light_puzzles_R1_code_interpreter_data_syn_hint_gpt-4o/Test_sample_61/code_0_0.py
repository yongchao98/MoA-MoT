# Initial grid
grid = [
    [1, 1],
    [1, 0]
]

# Button A controls these positions
positions = [(1, 1), (1, 0)]

# Toggle the lights controlled by Button A
for x, y in positions:
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Flatten the grid to a single list
result = [grid[i][j] for i in range(2) for j in range(2)]

# Print the result
print(result)