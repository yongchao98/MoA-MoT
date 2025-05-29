# Initial grid
grid = [
    [1, 0],
    [0, 0]
]

# Positions controlled by Button A (1-based index)
button_A_positions = [(2, 2), (2, 1)]

# Convert 1-based index to 0-based index and toggle the lights
for pos in button_A_positions:
    row, col = pos
    grid[row - 1][col - 1] = 1 - grid[row - 1][col - 1]

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(2) for j in range(2)]

# Print the final state
print(final_state)