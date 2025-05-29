# Initial grid
grid = [
    [1, 0],
    [0, 1]
]

# Positions controlled by Button A (1-based index)
button_a_positions = [(0, 0), (1, 0)]  # Converted to 0-based index

# Toggle the lights controlled by Button A
for pos in button_a_positions:
    x, y = pos
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(2) for j in range(2)]

# Print the final state
print(final_state)