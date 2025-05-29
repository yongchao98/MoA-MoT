# Initial grid configuration
grid = [
    [1, 0, 0],
    [1, 1, 1],
    [0, 1, 1]
]

# Button A controls these positions
button_A_positions = [(1, 2), (2, 0), (1, 1)]

# Button B controls these positions
button_B_positions = [(2, 1), (0, 0), (2, 0)]

# Function to toggle the lights at given positions
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Round 1: Press Button A
toggle_lights(grid, button_A_positions)

# Round 2: Press Button B
toggle_lights(grid, button_B_positions)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)