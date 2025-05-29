# Initial grid configuration
grid = [
    [0, 0, 1, 0, 1],
    [1, 0, 1, 1, 0],
    [0, 0, 1, 0, 0],
    [0, 1, 1, 0, 1],
    [0, 0, 1, 0, 1]
]

# Button control positions
button_A_positions = [(2, 5), (5, 5), (5, 2), (4, 3), (4, 4), (3, 1), (3, 2), (5, 3), (5, 1)]
button_B_positions = [(1, 4), (2, 4), (5, 1), (3, 5), (2, 3), (5, 2), (4, 5), (3, 4), (1, 1), (4, 4), (4, 3), (1, 5)]
button_C_positions = [(4, 2), (4, 5), (4, 1), (3, 5), (2, 1), (3, 4), (4, 3), (4, 4), (5, 2), (2, 4)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
toggle_lights(grid, button_C_positions)  # Round 1: Button C
toggle_lights(grid, button_A_positions)  # Round 2: Button A
toggle_lights(grid, button_A_positions)  # Round 3: Button A
toggle_lights(grid, button_A_positions)  # Round 4: Button A

# Flatten the grid to a single list
final_state = [cell for row in grid for cell in row]

# Print the final state
print(final_state)