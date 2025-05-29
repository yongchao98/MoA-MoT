# Initial grid state
grid = [
    [0, 1, 1, 0],
    [0, 1, 1, 1],
    [0, 1, 1, 0],
    [0, 1, 0, 1]
]

# Button control positions
button_A_positions = [(0, 1), (2, 0), (3, 0), (0, 2)]
button_B_positions = [(3, 1), (1, 1)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_A_positions)  # Round 1: Button A
toggle_lights(grid, button_B_positions)  # Round 2: Button B
toggle_lights(grid, button_A_positions)  # Round 3: Button A

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)