# Initial grid setup
grid = [
    [0, 1, 0, 1],
    [0, 1, 0, 0],
    [1, 0, 1, 1],
    [0, 0, 1, 0]
]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Button A positions (zero-based)
button_a_positions = [(0, 2), (1, 2), (3, 0), (3, 2)]

# Button B positions (zero-based)
button_b_positions = [(1, 1), (1, 2), (0, 1), (3, 0)]

# Round 1: Button A
toggle_lights(grid, button_a_positions)

# Round 2: Button B
toggle_lights(grid, button_b_positions)

# Round 3: Button A
toggle_lights(grid, button_a_positions)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)