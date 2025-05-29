# Initial grid state
grid = [
    [1, 1, 0, 0, 1],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1],
    [1, 1, 0, 0, 1],
    [1, 0, 0, 0, 1]
]

# Button control positions
button_A_positions = [(1, 3), (3, 2), (4, 3), (2, 4), (5, 1), (3, 3), (1, 4), (3, 1)]
button_B_positions = [(4, 3), (4, 2), (1, 4), (2, 4), (1, 5), (4, 5), (2, 5), (3, 1)]
button_C_positions = [(3, 1), (5, 2), (2, 4), (4, 3), (2, 5), (3, 4)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
toggle_lights(grid, button_A_positions)  # Round 1: Button A
toggle_lights(grid, button_A_positions)  # Round 2: Button A
toggle_lights(grid, button_A_positions)  # Round 3: Button A
toggle_lights(grid, button_B_positions)  # Round 4: Button B

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)