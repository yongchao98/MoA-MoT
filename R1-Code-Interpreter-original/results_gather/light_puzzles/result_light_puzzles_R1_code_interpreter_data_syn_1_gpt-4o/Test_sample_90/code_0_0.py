# Initial grid setup
grid = [
    [1, 0, 0, 0],
    [1, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 0, 0, 0]
]

# Button control mappings
button_A_positions = [(2, 3), (3, 4)]
button_B_positions = [(1, 2), (1, 4), (2, 1), (3, 1), (3, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Round 1: Button B
toggle_lights(grid, button_B_positions)

# Round 2: Button A
toggle_lights(grid, button_A_positions)

# Round 3: Button A
toggle_lights(grid, button_A_positions)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)