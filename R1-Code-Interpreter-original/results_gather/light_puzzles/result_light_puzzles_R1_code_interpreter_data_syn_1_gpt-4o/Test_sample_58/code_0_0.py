# Initial grid setup
grid = [
    [1, 0],
    [0, 0]
]

# Button A controls these positions
button_a_positions = [(0, 1), (1, 1)]  # zero-indexed positions

# Toggle the lights controlled by Button A
for x, y in button_a_positions:
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Flatten the grid to a single list for the final output
result = [grid[i][j] for i in range(2) for j in range(2)]
print(result)