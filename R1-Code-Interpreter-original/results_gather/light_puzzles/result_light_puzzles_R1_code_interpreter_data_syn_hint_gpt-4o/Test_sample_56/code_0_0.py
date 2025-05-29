# Initial grid
grid = [
    [0, 0],
    [1, 1]
]

# Positions controlled by Button A (1-based index)
button_a_positions = [(2, 1), (1, 1)]

# Convert to 0-based index and toggle the lights
for pos in button_a_positions:
    row, col = pos[0] - 1, pos[1] - 1
    grid[row][col] = 1 - grid[row][col]  # Toggle the light

# Flatten the grid and print the result
result = [grid[i][j] for i in range(2) for j in range(2)]
print(result)